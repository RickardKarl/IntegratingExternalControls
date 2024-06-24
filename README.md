# Robust integration of external control data in randomized trials

This repository is the official implementation of "Robust integration of external control data in randomized trials".

[Link to paper (available once public)]()

## Abstract

To be added

## Code instructions 

All experiments were run using R. 

- To see code used for the simulation study, see experiments/run_simulation_study.R. 
- Folders:
  - m_estimation: contains all estimators implemented using M-estimation, including our proposed randomization-aware and combined estimators.
  - dynamic_borrowing: contains alternative approaches based on a "dynamic borrowing" approach 
  - data: contains file used to generate data
  - experiments: contains files used to run simulations and generate plots.
  - output_folder: empty folder where results are saved after running experiment.
  - test: contains method to test the testable implication described in Sec. 4 of the paper.
  
  
  
#### Required dependencies
R 4.3.1 was used; below we list the necessary dependencies to run the code:
- MASS_7.3-60
- progress_1.2.3
- geex_1.1.1
- lmtest_0.9-40
- historicalborrow_1.0.4
- caret_6.0-94
- SelectiveIntegrative_1.0 (installed using devtools from Github repository: Gaochenyin/SelectiveIntegrative)
