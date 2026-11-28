## Gradient check: passed (known mismatches only)

3 programs, 56 points; 0 failed, 2 known mismatches, 0 skipped where log_prob is -inf.

| family | points | ok | skipped | failed | known | max rel. error vs central differences | note |
|---|---:|---:|---:|---:|---:|---:|---|
| cogmod_ddm | 22 | 20 | 0 | 0 | 2 | 6.4e-02 | wiener_lpdf() partials for w and sw are wrong at a wide st0; ndt far above the data at 1.2e-4 under the 1e-3 tolerance |
| cogmod_ddm4 | 16 | 16 | 0 | 0 | 0 | 5.9e-08 |  |
| cogmod_ddm5 | 18 | 18 | 0 | 0 | 0 | 1.7e-07 |  |

cogmod_ddm4: cogmod_ddm with sigmadrift = 0, sigmabias = 0, sigmandt = 0. cogmod_ddm5: cogmod_ddm with sigmabias = 0, sigmandt = 0.
