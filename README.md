# KFC-Procedure
KFC procedure is a three-step machine learning method aim at constructing predictions in both classification and regression problems.
It is available in the journal of Statistical Computation and Simulation at the following link: [https://dx.doi.org/10.1080/00949655.2021.1891539](https://dx.doi.org/10.1080/00949655.2021.1891539).

# Summary

This proceudure consists of three steps:
- K-step: K-means algorithm with M Bregman divergences. At the end of this step, we have M different partition structures of the input data.
- F-step: Fitting simple local models on the the K clusters obtained in the previous step for a given option of Bregman divergence. The collection of these K local models is called a candidate model. At the end of this step, we have M candidate models corresponding to M options of the Bregman divergences. For each Bregman divergence, we first assign the point to the closest center using the corresponding Bregman divergence, and the prediction is given by the corresponding local model on that cluster.
- C-step: Consensual aggregation method which will be used to combine all the predictions obtained from the M candidate models obtained in the previous step. The combining estimation methods (available in the other repository: https://github.com/hassothea/Kernel-based-Aggregation-Method_COBRA_MixCOBRA) are:
          - Kernel-based consensual regression aggregation method (regression)
          - Kernel-based consensual classification aggregation method (binary classification)
          - MixCOBRA (classification and regression).

The predictions given by individual candidate estimators and the final step of KFC procedure will be provided by the function.


