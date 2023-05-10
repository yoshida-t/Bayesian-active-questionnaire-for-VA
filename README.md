
<!-- README.md is generated from README.Rmd. Please edit that file -->

### introduction

This folder contains .R files for functions and experiments illustrated
in the paper “Bayesian active questionnaire design for cause-of-death
assignment using verbal autopsies” by T. Yoshida, T.S. Fan, T.
McCormick, Z. Wu & Z.R. Li. (2023) (The AHLI Conference on Health,
Inference, and Learning (CHIL) )

### File Structure

- rcodes
  - ‘functions.R’ : Contain all functions used for experiments
  - ‘0-simulate-data.R’ : Generate synthetic data (correctly specified
    model and misspecified model)
  - ‘1-experiment-nonstop.R’ : Conduct an experiment to compare accuracy
    over iterations with synthetic data
  - ‘2-experiment-stopping.R’ : Conduct an experiment to compare
    different stopping rules with synthetic data
  - ‘3-experiment-penalty.R’ : Conduct an experiment to compare
    different noise settings (‘h’) and different jumping penalization
    parameters (‘lambda’) with synthetic data
  - ‘4-experiment-phmrc-nonstop-CV.R’ : Conduct an experiment to compare
    accuracy over iterations with PHMRC data
  - ‘5-experiment-phmrc-stopping-CV.R’ : Conduct an experiment to
    compare different stopping rules with PHMRC data
  - ‘6-experiment-phmrc-stopping-process-CV.R’ : Conduct an experiment
    to compare different noise settings (‘h’) and different jumping
    penalization parameters (‘lambda’) with PHMRC data
  - cpp
    - ‘PWKL.cpp’: Run PWKL algorithm
  - ‘output’ : Store CV results
  - ‘table’ : Store tables for stopping rule experiments
  - ‘input_data’ : Store synthetic and PHMRC data
  - ‘fig’ : Store figures
  - ‘data’ : Dataset for PHMRC
