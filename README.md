# HCC prediction study 2023
This repository provides all the code used in Van Renne et al. for the manuscript: "A liver and serum IgA signature predicts hepatocellular 
carcinoma in chronic viral hepatitis patients"

This repository contains all the code to reproduce the gene signature, and figures of this manuscript.

All code was written for R. 

Some remarks for reproducing the data:

It is best to keep the folder structure as shown here.

The gene transcription data files are stored in NCBI GEO GSE237332 and processed data (RPM and raw read counts) can be retrieved from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237332
These gene expression matrices should be formatted as .gct files. More info on gct format on https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
The headers of the .gct files are provided in data_files/collapsed_data. The processed data from the training set and validation set on GSE237332 can be pasted below this header in order to work. 

The scRNA-seq data seurat object can be downloaded from Zenodo (DOI 10.5281/zenodo.10149894) as a .Rds file. It should be stored in data_files/scRNAseq_data 

In the output folder, all files will be stored that are created with this code and that are loaded in the different .R code modules. 
