# prophet
The simulation scripts are based on the method proposed in our paper "Personalization and Optimization of Decision Parameters via Heterogenous Causal Effects".

We share example scripts for conduct simulation analysis in exam- ining the proposed methods and stochastic optimization algorithms in the following Github link: https://github.com/tuye0305/prophet. The source code contains five main files:
(1) multipleRunSimulation.R - This code can generate multi- ple runs of simulation study flow (including generate simula- tion data, apply all proposed approaches to calculate optimal policy x∗ and evaluate the solution using simulated test datasets).
(2) generateSimulationDataUtils.R - This file have the util functions associated with generating simulation data with self-initiated causal DAG.
(3) cohortLevelProphetUtils.R - This file includes the func- tions to run cohort-level solution paired with stochastic op- timization (CT.ST).
(4) memberLevelProphetUtils.R - This file covers the util func- tions to enable the member-level estimations (both Causal Tree or Two-Model method with regression tree models) paired with deterministic optimization (CT.DT, TM.DT).
(5) StochasticOpSimulation.R - This code runs the stochastic optimization routine to identify the optimal parameter in each cohort. Running the code as is should generate the plot as shown in Figure 7.
We have also used the open source R libraries simcausal [17] to generate simulation datasets, and the Causal Tree library (https: //github.com/susanathey/causalTree) to identify the cohorts and estimate effects for each treatment j and metric k. Using these, the entire methodology discussed in this paper can be easily reproduced for any similar problem of interest for small scale problems.