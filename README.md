
IIIVmrMLM and mppQTL - - 

Please use the newest software (IIIVmrMLM_1.0.zip) uploaded on 2022-11-25.

If users want to identify QTN-by-environment interactions (QEIs), please use the IIIVmrMLM.QEI software. In this software, multi-environment joint analysis may be used to directly identify QTNs and QEIs (this function is also included in the IIIVmrMLM software), while single-environment analyses under random/fixed SNP-effect models, along with trait differences, regression intercept and coefficient, and environmental variation indicators (range, variance, standard deviation, and coefficient of variation) as phenotypes, are used to indirectly identify QEIs. The software is released on 2022-11-25.

In most existing methods and softwares of genome-wide association studies (GWAS) for detecting quantitative trait nucleotides (QTNs), QTN-by-environment interactions (QEIs), and QTN-by-QTN interactions (QQIs), only allele substitution effect and its interaction-related effects are detected and estimated, conditional on method-specific polygenic background control, leading to confounding in effect estimation and insufficient polygenic background control (Li et al., 2022a,b). To address these issues, we have recently established a compressed variance component mixed model method, 3VmrMLM, to estimate additive and dominant effects and their environmental and epistatic interaction effects conditional on fully controlling all possible polygenic backgrounds (Li et al., 2022a). To facilitate the implement of the 3VmrMLM method in GWAS, we developed IIIVmrMLM as the R and C++ softwares for 3VmrMLM (Li et al., 2022b).

Once the software ( version)  is installed, two commands below are available to detect QTNs:

library(IIIVmrMLM)
IIIVmrMLM(fileGen="D:/Users/Genotype",filePhe="D:/Users/Phenotype.csv",method="Single_env",trait=c(1:3),dir="D:/Users/")

Two commands below are available to detect QEIs:

library(IIIVmrMLM)

IIIVmrMLM(fileGen="D:/Users/Genotype", filePhe="D:/Users/Phenotype_multi_env.csv", method="Multi_env",trait=1:2,n.en=c(2,2),dir="D:/Users/")

Two commands below are available to detect QQIs:

library(IIIVmrMLM)

IIIVmrMLM(fileGen="D:/Users/Genotype",filePhe="D:/Users/Phenotype.csv",method="Epistasis",trait=c(1:3),dir="D:/Users/")

Some optional parameters in IIIVmrMLM() may be the default, such as fileKin, filePS, PopStrType, fileCov, SearchRadius, svpal, DrawPlot, Plotformat, and Chr_name_com. If users want to adjust these parameters, please refer to the software instruction for details.

If R version in your computer is R 4.2, please use this software IIIVmrMLM_1.0.zip. If R version in your computer is R 4.1, the previously-uploaded R software IIIVmrMLM_1.0.zip is also available. 

References

Li M, Zhang YW, Zhang ZC, Xiang Y, Liu MH, Zhou YH, Zuo JF, Zhang HQ, Chen Y, Zhang YM (2022a). A compressed variance component mixed model for detecting QTNs and QTN-by-environment and QTN-by-QTN interactions in genome-wide association studies. Molecular Plant 15(4): 630-650.

Li Mei, Zhang Ya-Wen, Xiang Yu, Liu Ming-Hui, Yuan-Ming Zhang. (2022b). IIIVmrMLM: The R and C++ tools associated with 3VmrMLM, a comprehensive GWAS method for dissecting quantitative traits. Molecular Plant 15(8): 1251-1253.


mppQTL --

This package implements a multi-locus method, namely mppQTL, to detect Quantitative trait loci (QTL) in multi-parent populations (MPPs), e.g., multiparent advanced generation intercross (MAGIC), nested association mapping (NAM), and random-open-parent association mapping (ROAM) populations. In addition, it allows users to customize design matrix (probabilities of genotypes) for identifying QTLs in any MPPs. The new method includes two steps. The first one is genome-wide scanning based on a single-locus linear mixed model, their P-values are obtained from likelihood ratio test, the peaks of negative logarithm P-value curve are further selected by group-lasso, and all the selected peaks are regarded as potential QTLs. In the second step, all the potential QTLs are placed on a multi-locus model, all the effects are estimated by expectation maximization empirical Bayes algorithm, and all the non-zero effects are further evaluated via likelihood ratio test for significant QTLs. This package provides an effective method of detecting small-effect QTLs in any MPPs. The methodological manuscript is in revision.
