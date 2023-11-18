# HCC prediction study 2023
This repository provides all the code used in Van Renne et al. for the manuscript: "A liver and serum IgA signature predicts hepatocellular 
carcinoma in chronic viral hepatitis patients"

This repository contains all the code to reproduce the gene signature, and figures of this manuscript.

All code was written for R. 

Some remarks for reproducing the data:

It is best to keep the folder structure as shown here.

The gene transcription data files are stored in NCBI GEO GSE237332 and processed data (RPM and raw read counts) can be retrieved from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237332
These gene expression matrices should be formatted as .gct files. More info on gct format on https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

The headers of the .gct files are provided in data_files/collapsed_data. The processed data from the training set and validation set on GSE237332 can be pasted in these files in order to work.  

The scRNA-seq data seurat object can be downloaded from Zenodo (DOI 10.5281/zenodo.10149894) as a .Rds file. It should be stored in data_files/scRNAseq_data 

In the output folder, all files will be stored that are created with this code and that are loaded in the different .R code modules. 


# order of code
1) calculate_signature.R
2) NTP_HCC_risk_prediction.R
3) volcanoplot_cox_scores.R
4) volcanoplot_diffexp.R
5) heatmaps_NTP.R
6) reproduce_scRNAseq.R (original data matrices can be requested to S. MacParland to reproduce. Output of this code can be retrieved from Zenodo as described above)
7) signatures_scRNAseq.R  
8) single_cell_deconvolution.R
9) IHC_staining_heatmap.R
10) serumIgA_bulkRNA_correlation.R
11) serum_cohort_IgA.R (all patients or patients with censored time>1 year can be set in this code)

    
