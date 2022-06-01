# KFC-Procedure
KFC procedure is a three-step machine learning method aim at constructing predictions in both classification and regression problems.
It is recently accepted and available in the journal of Statistical Computation and Simulation at the following link: https://dx.doi.org/10.1080/00949655.2021.1891539

# Abstract

Nowadays, many machine learning procedures are available on the shelve and may be used easily to calibrate predictive models on supervised data. However, when the input data consists of more than one unknown cluster, linked to different underlying predictive models, fitting a model is a more challenging task. We propose, in this paper, a three-step procedure to automatically solve this problem. The first step of the procedure aims at catching the clustering structure of the input data, which may be characterized by several statistical distributions. For each partition, the second step fits a specific predictive model based on the data in each cluster. The overall model is computed by a consensual aggregation of the models corresponding to the different partitions. A comparison of the performances on different simulated and real data assesses the excellent performance of our method in a large variety of prediction problems.

# Summary

This proceudure consists of three steps:
- K-step: K-means algorithm with M Bregman divergences. At the end of this step, we have M different partition structures of the input data.
- F-step: Fitting simple local models on the the K clusters obtained in the previous step for a given option of Bregman divergence. The collection of these K local models is called a candidate model. At the end of this step, we have M candidate models corresponding to M options of the Bregman divergences. For each Bregman divergence, we first assign the point to the closest center using the corresponding Bregman divergence, and the prediction is given by the corresponding local model on that cluster.
- C-step: Consensual aggregation method which will be used to combine all the predictions obtained from the M candidate models obtained in the previous step. The combining estimation methods (available in the other repository: https://github.com/hassothea/Kernel-based-Aggregation-Method_COBRA_MixCOBRA) are:
          - Kernel-based consensual regression aggregation method (regression)
          - Kernel-based consensual classification aggregation method (binary classification)
          - MixCOBRA (classification and regression).

The predictions given by individual candidate estimators and the final step of KFC procedure will be provided by the function.
