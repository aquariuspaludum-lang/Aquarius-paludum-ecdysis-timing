# README: Aquarius_paludum_ecdysis_timing.csv

---

## Overview

This dataset contains records of ecdysis timing in the water strider *Aquarius paludum* (Fabricius) (Hemiptera: Gerridae), collected to examine the effects of feeding timing and light regime on the daily timing of ecdysis. The data are associated with the following publication:

> Kishi M. (in press) Light-dependent effects of feeding timing on the daily timing of ecdysis in the water strider *Aquarius paludum* (Hemiptera: Gerridae). *[Journal — to be completed upon acceptance]*. https://doi.org/[to be added upon publication]

---

## Files included

| File | Description |
|---|---|
| `Aquarius_paludum_ecdysis_timing.csv` | Raw data: ecdysis events recorded three times daily |
| `01_statistical_analysis.R` | R script: log-linear model, chi-squared tests, synchrony index, Wilcoxon test |
| `02_gating_analysis.R` | R script: within-instar daily ecdysis distribution summary |

### How to run

Execute scripts in the following order. All files must be in the same working directory.

```r
source("01_statistical_analysis.R")   # all statistical tests
source("02_gating_analysis.R")        # within-instar shift summary
```

### Software requirements

| Software | Version |
|---|---|
| R | ≥ 4.0.0 |
| MASS | ≥ 7.3 |

---

## File information

| Item | Detail |
|---|---|
| Filename | `Aquarius_paludum_ecdysis_timing.csv` |
| Format | Comma-separated values (CSV), UTF-8 encoding |
| Number of rows | 240 (excluding header) |
| Number of columns | 7 |
| Total ecdysis events recorded | 3,307 individuals |

---

## Experimental overview

Adults of *A. paludum* were collected at the Kagami River, Kochi City (33°33'N, 133°30'E), Kochi Prefecture, Japan, in August 2004. Newly hatched first-instar larvae were assigned to one of four feeding groups under each of two light regimes (see below). All groups were established at the middle of the light phase. Ecdysis events from the 1st–2nd instar transition through the 5th–adult transition were recorded three times daily (just after light-on, middle of the light phase, and just before lights-off). The food provided was adults of the blowfly *Lucilia illustris* (Meigen); food and water were changed simultaneously with each feeding. Temperature was maintained at 20 ± 2°C throughout the experiment.

---

## Column descriptions

### `light_regime`
Light condition under which larvae were reared.

| Value | Condition |
|---|---|
| 1 | 15.5 h light : 8.5 h dark (15.5L : 8.5D) |
| 2 | 24 h continuous light (24L) |

---

### `feeding_group`
Daily feeding schedule assigned to each larval group. Groups under 15.5L : 8.5D are designated A–D; corresponding groups under 24L are designated A'–D'.

| Value | Group (15.5L:8.5D / 24L) | Feeding schedule |
|---|---|---|
| 1 | A / A' | Three times daily: just after light-on, middle of the light phase, and just before lights-off |
| 2 | B / B' | Once daily: just after light-on |
| 3 | C / C' | Once daily: middle of the light phase |
| 4 | D / D' | Once daily: just before lights-off |

Note: For group A / A' (feeding_group = 1), water was changed only at just after light-on; for all other groups, water was changed at the same time as feeding.

---

### `observation_sequence`
Sequential observation number assigned across the entire experiment, starting from the midday observation on the day of group establishment.

- Sequence 0: Day 1, middle of the light phase (group establishment day)
- Sequence 1: Day 1, just before lights-off
- Sequence ≥ 2: Day = (sequence − 2) ÷ 3 + 2 (integer division); time of day cycles as: just after light-on → middle of the light phase → just before lights-off

Note: The minimum value of `observation_sequence` in this dataset is 7, as no ecdysis events were observed during the initial observations (sequences 0–6).

---

### `observation_day`
Calendar day of the observation, derived from `observation_sequence` according to the rule described above. Day 1 corresponds to the day of group establishment.

---

### `observation_time`
Time of day at which the ecdysis event was recorded.

| Value | Observation time |
|---|---|
| 1 | Just after light-on (morning; approximately ZT 0) |
| 2 | Middle of the light phase (midday; approximately ZT 7.75) |
| 3 | Just before lights-off (evening; approximately ZT 15.5) |

Note: ZT (Zeitgeber time) is defined with ZT 0 = lights-on under 15.5L:8.5D. Under 24L (continuous light), observations were conducted at the same clock times as the 15.5L:8.5D groups; ZT is not applicable under constant light conditions. Within each `observation_sequence`, the value of `observation_time` is unique (one-to-one correspondence).

---

### `instar`
Larval instar *entered* at the time of the recorded ecdysis event.

| Value | Stage |
|---|---|
| 2 | 2nd larval instar (ecdysis from 1st to 2nd instar) |
| 3 | 3rd larval instar (ecdysis from 2nd to 3rd instar) |
| 4 | 4th larval instar (ecdysis from 3rd to 4th instar) |
| 5 | 5th larval instar (ecdysis from 4th to 5th instar) |
| 6 | Adult emergence (ecdysis from 5th instar to adult) |

---

### `n_ecdysis`
Number of individuals recorded undergoing ecdysis at the given observation time, for the specified combination of `light_regime`, `feeding_group`, `instar`, and `observation_sequence`.

---

## Summary statistics

| light_regime | feeding_group | n_ecdysis (total) |
|---|---|---|
| 1 (15.5L:8.5D) | 1 (A) | 424 |
| 1 (15.5L:8.5D) | 2 (B) | 409 |
| 1 (15.5L:8.5D) | 3 (C) | 374 |
| 1 (15.5L:8.5D) | 4 (D) | 406 |
| 2 (24L) | 1 (A') | 429 |
| 2 (24L) | 2 (B') | 424 |
| 2 (24L) | 3 (C') | 432 |
| 2 (24L) | 4 (D') | 409 |
| **Total** | | **3,307** |

Ecdysis events are recorded across five instar transitions (2nd through 6th/adult) and 240 combinations of light regime, feeding group, instar, and observation sequence.

---

## Notes

- Food ration per feeding: 1 fly per 5 first- or second-instar larvae; 1 fly per 3 third- or fourth-instar larvae; 1 fly per 2 fifth-instar larvae (Harada, 1992).
- Rows in the dataset represent unique combinations of `light_regime`, `feeding_group`, `instar`, and `observation_sequence` for which at least one ecdysis event was recorded (`n_ecdysis` ≥ 1). Combinations with zero ecdysis events are not included.

---

## Contact

Manabu Kishi  
Persimmon and Peach Laboratory, Wakayama Fruit Tree Experiment Station  
Kinokawa, Wakayama 649-6531, Japan  
E-mail: kishi@hotmail.co.jp  
ORCID: 0009-0001-8955-5189

---

## License

This dataset is released under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.  
You are free to share and adapt the material for any purpose, provided appropriate credit is given.  
https://creativecommons.org/licenses/by/4.0/

---

## Citation

If you use this dataset, please cite the associated paper:

> Kishi M. (in press) Light-dependent effects of feeding timing on the daily timing of ecdysis in the water strider *Aquarius paludum* (Hemiptera: Gerridae). *[Journal — to be completed]*. https://doi.org/[to be added]

Data archived at: GitHub. https://github.com/aquariuspaludum-lang/Aquarius-paludum-ecdysis-timing (DOI to be added upon publication)

---

## References

Harada, T. (1992). The oviposition process in two direct breeding generations in a water strider, *Aquarius paludum* (Fabricius). *Journal of Insect Physiology*, 38(9), 687–692. https://doi.org/10.1016/0022-1910(92)90050-N
