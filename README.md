
IIIVmrMLM - - 

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

References

Li M, Zhang YW, Zhang ZC, Xiang Y, Liu MH, Zhou YH, Zuo JF, Zhang HQ, Chen Y, Zhang YM (2022a). A compressed variance component mixed model for detecting QTNs and QTN-by-environment and QTN-by-QTN interactions in genome-wide association studies. Molecular Plant 15(4): 630-650.

Li Mei, Zhang Ya-Wen, Xiang Yu, Liu Ming-Hui, Yuan-Ming Zhang. (2022b). IIIVmrMLM: the R and C++ tools associated with 3VmrMLM, a comprehensive GWAS method for dissecting quantitative traits. Molecular Plant, online (2022-06-08), DOI: 10.1016/j.molp.2022.06.002.
