# Robust integration of external control data in randomized trials

[Link to paper](https://arxiv.org/abs/2406.17971)

## Abstract

One approach for increasing the efficiency of randomized trials is the use of "external controls" -- individuals who received the control treatment in the trial during routine practice or in prior experimental studies. Existing external control methods, however, can have substantial bias if the populations underlying the trial and the external control data are not exchangeable. Here, we characterize a randomization-aware class of treatment effect estimators in the population underlying the trial that remain consistent and asymptotically normal when using external control data, even when exchangeability does not hold. We consider two members of this class of estimators: the well-known augmented inverse probability weighting trial-only estimator, which is the efficient estimator when only trial data are used; and a more efficient member of the class when exchangeability holds and external control data are available, which we refer to as the optimized randomization-aware estimator. To achieve robust integration of external control data in trial analyses, we then propose a combined estimator based on the efficient trial-only estimator and the optimized randomization-aware estimator. We show that the combined estimator is consistent and no less efficient than the most efficient of the two component estimators, whether the exchangeability assumption holds or not. We examine the estimators' performance in simulations and we illustrate their use with data from two trials of paliperidone extended-release for schizophrenia.

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
