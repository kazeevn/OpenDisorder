## Synthetic cases

| case                                        | expected   | MatterGen   | opendisorder           |
|:--------------------------------------------|:-----------|:------------|:-----------------------|
| LiFeO2 ordering of (Li,Fe)O                 | True       | True        | disordered_realization |
| ordered vacancies LiO2 of Li0.5O            | True       | False       | disordered_realization |
| Li3FeO4 (off-stoichiometry)                 | False      | False       | disordered_family      |
| vacancies filled: LiO                       | False      | False       | disordered_family      |
| missing O sublattice: LiFe                  | False      | False       | none                   |
| half-empty cation lattice: LiO2 vs (Li,Fe)O | False      | False       | none                   |
| foreign element Co                          | False      | False       | none                   |
| different framework (zinc blende)           | False      | False       | none                   |

## Literature cases

| opendisorder            |   match |   none |
|:------------------------|--------:|-------:|
| disordered_family       |       0 |      3 |
| disordered_realization  |      13 |      0 |
| none                    |       0 |     43 |
| pristine                |       1 |      0 |
| substitutional_ordering |       0 |     12 |
