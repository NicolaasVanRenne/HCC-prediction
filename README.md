# HCC prediction study 
This repository provides all the code used in Van Renne et al. for the manuscript: "A liver and serum IgA signature predicts hepatocellular 
carcinoma in chronic viral hepatitis patients"

This repository contains all the code to reproduce the gene signature, and figures of this manuscript.


Some remarks for reproducing the data:

It is best to keep the folder structure as shown here.

The gene transcription data files are stored in NCBI GEO GSE237332 and processed data (RPM and raw read counts) can be retrieved from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237332
These gene expression matrices should be formatted as .gct files. More info on gct format on https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

The headers of the .gct files are provided in data_files/collapsed_data. The processed data from the training set and validation set on GSE237332 can be pasted in these files in order to work.  

The scRNA-seq data seurat object can be downloaded from Zenodo (DOI 10.5281/zenodo.10149894) as a .Rds file. It should be stored in data_files/scRNAseq_data 

In the output folder, all files will be stored that are created with this code and that are loaded in the different .R code modules. 


## order of code
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

    




Installation
------------
All code was written for R, and can be easily run in the R environment. 

## Dependencies 
* R version >= 4.0.3. 
* R packages: survival, data.table, survminer, patchwork, dplyr, ggrepel, ggplot2, circlize, ComplexHeatmap, Seurat, harmony, ggsci, RColorBrewer, MuSiC, Biobase, ggpubr, lubridate, pROC

The code was run in R x64 4.0.3 installed on operating system Windows 10

## Issues
All feedback, bug reports (if any) and suggestions are welcome

## details of the R environment used by calling sessionInfo()

``` r
sessionInfo()

R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)
	
Matrix products: default

locale:
[1] LC_COLLATE=Dutch_Belgium.1252  LC_CTYPE=Dutch_Belgium.1252    LC_MONETARY=Dutch_Belgium.1252 LC_NUMERIC=C                   LC_TIME=Dutch_Belgium.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     
	
other attached packages:
[1] ggplot2_3.4.0      ggsci_2.9          harmony_0.1.0      Rcpp_1.0.7         SeuratObject_4.0.4 Seurat_4.0.2      
	
loaded via a namespace (and not attached):
	  [1] nlme_3.1-149          matrixStats_0.58.0    spatstat.sparse_2.0-0 RcppAnnoy_0.0.18      RColorBrewer_1.1-2    httr_1.4.2            sctransform_0.3.2     tools_4.0.3          
	  [9] utf8_1.2.1            R6_2.5.0              irlba_2.3.3           rpart_4.1-15          KernSmooth_2.23-17    uwot_0.1.10           mgcv_1.8-33           DBI_1.1.1            
	 [17] lazyeval_0.2.2        colorspace_2.0-0      withr_2.5.0           tidyselect_1.1.0      gridExtra_2.3         compiler_4.0.3        cli_3.4.1             plotly_4.9.3         
	 [25] labeling_0.4.2        scales_1.2.1          lmtest_0.9-38         spatstat.data_2.1-0   ggridges_0.5.3        pbapply_1.4-3         goftest_1.2-2         stringr_1.4.0        
	 [33] digest_0.6.27         spatstat.utils_2.3-1  pkgconfig_2.0.3       htmltools_0.5.2       parallelly_1.32.0     fastmap_1.1.0         htmlwidgets_1.5.3     rlang_1.0.6          
	 [41] shiny_1.6.0           farver_2.1.0          generics_0.1.0        zoo_1.8-9             jsonlite_1.7.2        ica_1.0-2             dplyr_1.0.5           magrittr_2.0.1       
	 [49] patchwork_1.1.1       Matrix_1.3-4          munsell_0.5.0         fansi_0.4.2           abind_1.4-5           reticulate_1.18       lifecycle_1.0.3       stringi_1.5.3        
	 [57] MASS_7.3-53           Rtsne_0.15            plyr_1.8.6            grid_4.0.3            parallel_4.0.3        listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1        
	 [65] crayon_1.4.1          miniUI_0.1.1.1        deldir_1.0-6          lattice_0.20-41       cowplot_1.1.1         splines_4.0.3         tensor_1.5            pillar_1.6.0         
	 [73] igraph_1.2.6          spatstat.geom_2.4-0   future.apply_1.7.0    reshape2_1.4.4        codetools_0.2-16      leiden_0.3.7          glue_1.4.2            data.table_1.14.0    
	 [81] png_0.1-7             vctrs_0.5.1           httpuv_1.5.5          gtable_0.3.0          RANN_2.6.1            purrr_0.3.4           spatstat.core_2.0-0   polyclip_1.10-0      
	 [89] tidyr_1.1.3           scattermore_0.7       future_1.26.1         assertthat_0.2.1      mime_0.10             xtable_1.8-4          RSpectra_0.16-0       later_1.1.0.1        
	 [97] survival_3.2-7        viridisLite_0.3.0     tibble_3.1.0          cluster_2.1.0         globals_0.15.1        fitdistrplus_1.1-3    ellipsis_0.3.2        ROCR_1.0-11          
```


