Efficiency Model Requirements
====

API
----
Implementation in python => expose a python API:

```
efficiency(k_p, pi1_p, pi2_p, pi3_p, phase_difference, decay_time, kaon_charge, year, magnet)
```

### Params
| Parameter   | Requirements |    Meaning       |
| ---------   | ------------ | --------------   |
|    `k_p`    | length 4, floats     | kaon 4-momentum  |
|    `pi1_p`  | length 4, floats     | pi1 4-momentum. <br> pi1 defined by  (Kpi1) < M(Kpi2) |
|    `pi2_p`  | length 4, floats     | pi2 4-momentum. <br> pi2 defined by (Kpi1) < M(Kpi2)  |
|    `pi3_p`  | length 4, floats     | pi3 4-momentum. <br> pi3 has same charge as the kaon  |
| `decay_time` | float | Proper decay time of the D |
| `kaon_charge` | +-1 | Charge of the Kaon |
| `year` | Int in {2011..2018} | Data taking year |
| `magnet` | String in {"MagUp", "MagDown"} | Magnetisation direction |


### Returns
Estimated detector efficiency for the provided event, as a float.

### Implementation
 - Prerequisite: trained BDT for each data taking year/magnetisation/phsp bin, each pickled
 - Look up the right BDT based on k charge, magnetisation, year, phsp bin
 - Use it to return a weight
 - BDT should be trained based on decay time and on phsp
