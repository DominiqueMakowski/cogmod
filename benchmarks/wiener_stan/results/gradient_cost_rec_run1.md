## Gradient cost: no regression (worst ratio 0.77)

Base 0.3.5 @ 6463109 vs PR 0.3.5 @ 6463109. Milliseconds per `grad_log_prob()` on 5000 simulated trials, median of 21 alternating blocks of 25; ratio = PR / base, flagged above 1.10, failing above 1.30.

| family | base ms | PR ms | ratio | |
|---|---:|---:|---:|---|
| cogmod_ddm4 | 12.165 | 9.358 | 0.77 | faster |
| cogmod_ddm5 | 35.348 | 10.626 | 0.30 | faster |

cogmod_ddm4: cogmod_ddm with sigmadrift = 0, sigmabias = 0, sigmandt = 0.
