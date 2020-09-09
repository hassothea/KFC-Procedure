# KFC-Procedure
KFC procedure is a three-step machine learning procedure used to build a predictive model in both classification and regression problems.
Read more: https://arxiv.org/abs/1909.09370

# Summary
This proceudure consists of three steps:
- K-step: K-mean algorithm with $M$ Bregman divergences. Each of them would lead to $M$ different partition structures of the input data.
- F-step: Fitting simple local models (linear regression for regression and logistic regression for binary classification problem) on the the $K$ clusters obtained in the previous step for a given option of Bregman divergence. The collection of these $K$ local models is called a candidate model. At the end of this step, we have $M$ candidate models corresponding to $M$ options of the Bregman divergences.
- C-step: Consensual aggregation method which will be used to combine all the $M$ candidate models obtained in the previous step. The combining estimation method are:
          - Kernel-based COBAR (regression)
          - Kernel-based consensual classification aggregation method (binary classification)
          - MixCOBRA (classification and regression).
