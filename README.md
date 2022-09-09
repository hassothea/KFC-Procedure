# KFC-Procedure

KFC procedure is a three-step machine learning method aim at constructing predictions in both classification and regression problems.
It is available in the journal of Statistical Computation and Simulation at the following link: [https://dx.doi.org/10.1080/00949655.2021.1891539](https://dx.doi.org/10.1080/00949655.2021.1891539).

--- 

## Summary

This proceudure consists of three steps:
- K-step: K-means algorithm with $M$ Bregman divergences (see, for example, [Banerjee (2005)](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://www.jmlr.org/papers/volume6/banerjee05b/banerjee05b.pdf)). Depending on the theoretical property of Bregman divergences, different Bregman divergences may cluster the input data into different structures. At the end of this step, we have M different partition structures of the input data.
- F-step: Fitting simple local models on the the K clusters obtained in the previous step for a given option of Bregman divergence. The collection of these K local models is called a candidate model. At the end of this step, we have M candidate models corresponding to M options of the Bregman divergences. For each Bregman divergence, we first assign the point to the closest center using the corresponding Bregman divergence, and the prediction is given by the corresponding local model on that cluster.
- C-step: Consensual aggregation method which will be used to combine all the predictions obtained from the M candidate models obtained in the previous step. The combining estimation methods (available in [AggregationMethods](https://github.com/hassothea/AggregationMethods) repository) are:
   - [`KernelAggReg`](https://hassothea.github.io/files/CodesPhD/KernelAggReg.html) : A kernel-based consensual regression aggregation method (see [Has (2021)](https://hal.archives-ouvertes.fr/hal-02884333v5)).
   - [`MixCobraReg`](https://hassothea.github.io/files/CodesPhD/MixCobraReg.html) : Aggregation using input-output trade-off (see [Fischer and Mougeot (2019)](https://www.sciencedirect.com/science/article/pii/S0378375818302349)).
   - [`KernelAggClass`](https://hassothea.github.io/files/CodesPhD/KernelAggClass.html) : A kernel-based combined classification rule (see [Mojirsheibani (2021)](https://www.sciencedirect.com/science/article/pii/S0167715200000249)).
   - [`MixCobraClass`](https://hassothea.github.io/files/CodesPhD/MixCobraClass.html) : Aggregation using input-output trade-off (see [Fischer and Mougeot (2019)](https://www.sciencedirect.com/science/article/pii/S0378375818302349)).

The predictions given by individual candidate estimators and the final step of KFC procedure will be provided by the function.

![The summary of KFC procedure.](kfc.png)

--- 

## &#128270; How to download & run the codes?

To run the codes, you can <span style="color: #097BC1">`clone`</span> the repository directly or simply load the <span style="color: #097BC1">`R script`</span> source files from this repository using [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package in `Rstudio` as follows:

1. Install [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package using command: 

    `install.packages("devtools")`

2. Loading the source codes from `GitHub` repository using `source_url` function by: 

    `devtools::source_url("https://raw.githubusercontent.com/hassothea/KFC-Procedure/master/file.R")`

where `file.R` is the file name contained in this repository which you want to import into your `Rstudio`.

---

## &#128214; Documentation

The documentation and explanation of the methods are available on my webpage as listed below:

- `KFCReg` : see [KFCReg documentation](https://hassothea.github.io/files/CodesPhD/KFCReg.html).

----
