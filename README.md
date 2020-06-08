# Dissertation
**Title:** Use of aggregated covariates in propensity score analysis of clustered data  
**Author:** Kyle Nickodem

Consider the increased presence of school resource officers (SROs) in American schools over the last 20 years. SROs are oft-armed law enforcement agents tasked with deterring misbehavior and increasing the students’ and staffs’ sense of safety in order to create a school environment conducive for learning.

*How might we determine whether the presence of an SRO has an effect on students' perceptions of school safety and support?*

Ideally, we conduct a Randomized Controlled Trial (RCT), but three issues complicate such a study:
1. Ethical and logistical feasibility of randomly assigning SROs to each school (i.e., cluster)
2. Accounting for the clustering of student outcomes within schools
3. Availability of school-level information

Propensity score methods address issue #1; multilevel models address issue #2; student information aggregated by school addresses issue #3. Great! So, if I'm an applied researcher:

*How do I conduct a propensity score analysis using multilevel models and aggregated covariates to evaluate the impact of a cluster-level treatment on subject-level outcomes?*

I conducted two simulation studies to explore study design choices on treatment effect estimation. I then utilized student responses from a statewide data collection to investigate the association between SRO presence in schools and student outcomes.

## Simulation 1

Various procedures for creating aggregated covariates have been used to evaluate their use in multilevel models, analysis of survey data, and in measurement models. However, it is unclear whether these procedures are suitable for the context of propensity score methods.

Questions:
1. What sample characteristics are generated from four procedures for creating aggregated covariates?
2. Are the generated samples suitable for a propensity score analysis of clustered data?

Files:
- Simulation Script - Sim1GenAgg.R
- Results Data - Sim1_ThousandRepsUpdate.R
- Analysis Script - Sim1Analysis.R

**Follow-up simulation**

Question: How many clusters are needed to minimize model non-convergence when estimating propensity scores?
- File - Sim1ConvergenceTest.R

## Simulation 2

In an RCT, treatment and control groups are created via random assignment. Random assignment can occur by clusters or by subjects. A propensity score analysis is defined by the lack of random assignment. Consequently, it is sometimes unclear in propensity score studies if treatment exposure should be appraised by clusters or by subjects. For instance. should the propensity score be defined as "the probability of a school having an SRO" (i.e., appraisal by clusters) or as "the probability of a student attending a school with an SRO" (i.e., appriasal by subjects)?

Additionally, a true confounder is a variable associated with both treatment exposure and the outcome. When a true confounder is missing from a propensity score analysis, treatment effect estimation can be biased. When a true confounder at the cluster level is missing, can it be replaced by a correlated aggregated covariate? For instance, if school-level socioeconomic status is a true confounder, is students' free/reduced price lunch status aggregated by school a viable replacement? 

Questions:
1. What impact, if any, does estimating the propensity score by cluster or by subjects have on covariate balance and treatment effect estimation?
2. What impact, if any, does replacing  true cluster-level confounders with correlated aggregated covariates have on covariate balance and treatment effect estimation?

Files:
- Simulation Script - Sim2ClusterPS.R
- Results Data - Sim2_ThousandRep_Pooled.R
- Analysis Script - Sim2Analysis.R

## Empirical Study

How well do the lessons learned from Simulation 2 translate to real world data? The practical problem addressed is the association between SRO presence in a school and students' preceptions of school safety and support.

Questions:
1. What is the association between SRO presence and students' perceptions of school safety and support?
2. How is the association impacted, if at all, by appriasing SRO presence at the school level or at the student level in the propensity score analysis?

Files:
- Data Cleaning Script - Applied_DataPrep.R
- Propensity Score Matching Script - Applied_PSM.R
- Outcome Analysis Script - Applied_HLM.R

**Sensitivity Analysis**

Does removing aggregated covariates from the analysis change the estimated association?
- File - Applied_AggSensitivity.R
