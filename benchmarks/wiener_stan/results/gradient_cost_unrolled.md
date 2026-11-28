## Gradient cost: no regression (worst ratio 0.80)

Base 0.3.5 @ 6463109 vs PR 0.3.5 @ 6463109. Milliseconds per `grad_log_prob()` on 5000 simulated trials, median of 21 alternating blocks of 25; ratio = PR / base, flagged above 1.10, failing above 1.30.

| family | base ms | PR ms | ratio | |
|---|---:|---:|---:|---|
| cogmod_ddm4 | 9.761 | 7.789 | 0.80 | faster |
| cogmod_ddm5 | 36.966 | 12.525 | 0.34 | faster |

cogmod_ddm4: cogmod_ddm with sigmadrift = 0, sigmabias = 0, sigmandt = 0.
