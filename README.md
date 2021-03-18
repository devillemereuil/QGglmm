# QGglmm

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/QGglmm)](https://cran.r-project.org/package=QGglmm)

## NEWS 

Due to its removal from CRAN, QGglmm dropped R2Cuba as a dependency to solve multivariate integrals. It is now using the package [cubature](https://CRAN.R-project.org/package=cubature). By taking advantage of the ["vectorised"](https://CRAN.R-project.org/package=cubature/vignettes/cubature.html) version of the algorithm, the multivariate computations of QGglmm (QGmvparams, QGvcov, QGmvmean, QGmvpsi, QGmvicc, QGmvpred) are considerably faster. Most functions are 10x-50x faster, but especially QGmvicc is 100x-500x faster. A comparison between the old and new version of the example of the man page of QGmvicc showed a decreased in computation from 25 minutes to... 4 seconds!

## What is this package?

QGglmm computes various quantitative genetics parameters on the observed data scale from latent parameters estimated using a Generalised Linear Mixed Model (GLMM) estimates. Especially, it yields the phenotypic mean, phenotypic variance and additive genetic variance on the observed data scale.

More information can be found in this [article](http://dx.doi.org/10.1534/genetics.115.186536) and on [CRAN](https://CRAN.R-project.org/package=QGglmm).

## How to install this package

### Using CRAN
* Simply use `install.packages("QGglmm")` as for any package.

### From this GitHub

* Install the packages on which QGglmm depends: cubature and mvtnorm.
    `install.packages(c("cubature","mvtnorm"))`
* Go to the [release page](https://github.com/devillemereuil/QGglmm/releases) and download the latest release.
* In a terminal, go to the folder where the release was downloaded and enter the following line:  
    `R CMD INSTALL QGglmm-xx.tar.gz` where `xx` is the version number.
* Alternatively, you can use the graphical tools of R-GUI or RStudio to manually install the package after download. For RStudio, this can be done using "Install Packages..." in the Tools menu, choosing "Install from: Package Archive File".

## Submit feedback

If you encounter any bug or usability issue, or if you have some suggestions or feature request, please use the [issue tracker](https://github.com/devillemereuil/QGglmm/issues). Thank you!
