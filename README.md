# Estimation of natural selection and allele age from time series allele frequency data using a novel likelihood-based approach
The source code is implemented for a likelihood-based method for inferring natural selection and allele age from time series data of allele frequencies, and the article has been published in Genetics, available at https://doi.org/10.1534/genetics.120.303400.

[Code v1.0](https://github.com/zhangyi-he/WFM-1L-DiffusApprox-KBE/tree/master/Code%20v1.0) includes the source code implemented for the case of constant demographies.

[Code v1.1](https://github.com/zhangyi-he/WFM-1L-DiffusApprox-KBE/tree/master/Code%20v1.1) includes the source code implemented for the case of non-constant demographies.

There have been some inquiries regarding how to run our original R code, which is basically not understandable. We therefore release a web app (a javascript/C++ implementation of the method outlined in https://doi.org/10.1534/genetics.120.303400) for the public (see https://fy-bris.github.io/saa-mle). Our app has been tested in Google Chrome version 101.0.4951.64, but the performance on any other browser is not guaranteed. A separate command-line program will be available soon.
 
Please let us know if you have any feedback. We will correct any bugs in this program in a month or two, and release an updated version.
