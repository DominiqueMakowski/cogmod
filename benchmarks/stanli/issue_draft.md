**Title:** Parameter-dependent branches in brms custom families: `log_mix()` gap, comparisons as values, and the cost of the region path

Hi, and thanks for stanli. I maintain [cogmod](https://github.com/DominiqueMakowski/cogmod), an R package of brms custom families for cognitive models (reaction-time distributions, race and diffusion models). I tried the light route from paul-buerkner/brms#1911 (`brm(empty = TRUE)` → `stanli_model()` → `sample_model()` → `as_stanfit()`) on one of our models: a Log-Normal Race with an outlier mixture, 4620 trials. The short version: once rewritten without branches, it runs **3-4x faster than CmdStan** with an identical posterior, which is great. Getting there took some rewriting that I suspect other custom-family authors will hit too. Everything below is reproduced with stanli 0.18.1 alone in [`repro_issue.R`](https://github.com/DominiqueMakowski/cogmod/blob/dev/benchmarks/stanli/repro_issue.R); the full harness is [here](https://github.com/DominiqueMakowski/cogmod/tree/dev/benchmarks/stanli).

### 1. `log_mix()` inside a parameter-dependent branch does not compile

```stan
functions {
  real f(real y, real mu) {
    if (y - mu <= 0) return log(0.1);
    return log_mix(0.1, -1.0, normal_lpdf(y | mu, 1));
  }
}
```

gives `stanli compile: runtime-control region: function log_mix`. Without the `if` it compiles, and so does the same `if` with `log_sum_exp(log(0.1) - 1.0, log1m(0.1) + normal_lpdf(y | mu, 1))`. All our densities end in `log_mix()` behind an early return like this, so this one call is what stops the brms-generated programs from compiling at all. `erfc`, `log1p`, `log1m_exp`, `lognormal_lpdf()` and calls to user functions all compiled fine in the same position.

### 2. A comparison used as a value does not compile

`return log(y > mu);` gives `unsupported function Greater__`, with or without a branch around it. `y > mu ? a : b`, `if (y > mu)` and `step(y - mu)` all compile.

### 3. A parameter-dependent branch is expensive, even when never taken

The same function with and without a far-tail switch, N = 5000. Every `y` is near `mu`, so the branch is never taken: same `lp`, same gradient.

```stan
// no_branch
return log(0.5 * erfc(-(y - mu) * 0.7071067811865476));
// branch
if (y - mu < -25) return -0.5 * square(y - mu) - log(mu - y) - 0.9189385332046727;
return log(0.5 * erfc(-(y - mu) * 0.7071067811865476));
```

`log_prob_grad()`: 750 µs without the branch, 2000 µs with it (2.7x). On the real model it compounds: the brms program with only `log_mix()` written out (branches kept) costs **27 ms** per gradient against **4.2 ms** in CmdStan 2.38 (Windows, i7-1265U).

We cannot drop the branches; they are the numerics. One is a tail switch so `log Phi` keeps a finite gradient past `x = -38`; without it, one slow trial turns the whole gradient into NaN. Another is the `t <= ndt` early return. So we rewrote every branch as a select:

```stan
real xl = fmin(x, -25);  real xh = fmax(x, -25);   // each arm clamped into its own domain
real w = step(-25 - x);
return w * lo(xl) + (1 - w) * hi(xh);              // the unpicked arm is finite, weight 0, zero adjoint
```

and the early return as a mask inside the one `log_mix()`: `log_mix(p, lp_out, lp_dec + (step(Y - ndt) - 1) * 1e300)`. That version matches CmdStan to 10 significant digits in `lp` and ~1e-15 in the gradient, including at `x = -65`, and costs **1.35 ms**: 3.1x faster than CmdStan, 3.7x per chain in a 4-chain fit.

### 4. Two performance cliffs we could only see on the full program

On the 4620-trial LNR, with no parameter-dependent branch anywhere:
- `log_sum_exp(log(p) + a, log1m(p) + b)` in place of `log_mix(p, a, b)`: 4.7 ms against 0.9.
- Two early returns on **data only** (`if (Y <= 0) return negative_infinity();`), which I expected to be folded at build time: 4.8 ms against 1.1.
- The `t <= ndt` switch written as a blend, `w * log_mix(...) + (1 - w) * (...)`, rather than a mask inside `log_mix()`: 5.5 ms against 0.9.

In a minimal program the first two cost 1.2x and nothing, so they need the surrounding program. [`bisect.R`](https://github.com/DominiqueMakowski/cogmod/blob/dev/benchmarks/stanli/bisect.R) builds all the variants from the emitted brms program, if useful.

### Suggestions, in rough order of how much they would help us

1. If-conversion of small, pure, two-armed branches into the select above, automatically. Code written the ordinary way would then stay on the graph, and nobody would need to rewrite it by hand.
2. `log_mix()` (and generally the `log_*` combinators) inside runtime-control regions.
3. Folding early returns on data, wherever the cost in 4 comes from.
4. Comparisons as real-valued expressions (`Greater__` and friends), like `step()`.
5. Minor, R package: `sample_model(init = )` with a list of per-chain init lists, as cmdstanr takes. Ours errors with "'list' object cannot be coerced to type 'double'".

Happy to test anything against the harness.
