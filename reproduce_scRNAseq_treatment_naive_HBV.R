################################################
# Narmada data set with annotated plasma cells
################################################

#set wd
	setwd("/your_data_path")

#load libraries
	library(Seurat)
	library(ggplot2)
	library(RColorBrewer)
	library(ggsci)
	#library(patchwork)

# set save.plots to save output files
	save.plots=FALSE

# Description of HBV single cell data set: Van Renne et al. 2025 PLoS One [PMID:39932940] Detection of hepatitis B virus mRNA from single cell RNA sequencing data without prior knowledge
# download seurat object from Zenodo "Detection of HBV mRNA in a liver scRNA-seq data set": https://zenodo.org/records/13643041
	SO.integrated <- readRDS(file="/mnt/bctl/nvan0112/data/scRNAseq_data/VanRenne2025/HBV_scRNAseq.Rds")

#Subset and annotate plasma cells
	SO.integrated -> SO.backup
	SO.integrated <- subset(SO.integrated, idents="plasmocyte")
	
#rename plasma cells according to heavy chain expression
	Idents(SO.integrated, cells = WhichCells(SO.integrated, expression = `IGHG1` > 4))  <- "IgG plasma cell"
	Idents(SO.integrated, cells = WhichCells(SO.integrated,  expression = `IGHG2` > 4)) <- "IgG plasma cell"
	Idents(SO.integrated, cells = WhichCells(SO.integrated,  expression = `IGHG3` > 4)) <- "IgG plasma cell"
	Idents(SO.integrated, cells = WhichCells(SO.integrated,  expression = `IGHG4` > 4)) <- "IgG plasma cell"
	Idents(SO.integrated, cells = WhichCells(SO.integrated,  expression = `IGHM` > 4))  <- "IgM plasma cell"
	Idents(SO.integrated, cells = WhichCells(SO.integrated,  expression = `IGHA1` > 4)) <- "IgA plasma cell"
	Idents(SO.integrated, cells = WhichCells(SO.integrated,  expression = `IGHA2` > 4)) <- "IgA plasma cell"
	
	#investigate remaining cells not annotated yet
		VlnPlot(SO.integrated, "nCount_RNA") 

		dotplot.features <- c("MZB1","TNFRSF17","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","IGHM","IGHD","BANK1","MS4A1","TNFRSF13C")
		p=DotPlot(SO.integrated, cols = "RdYlBu",features = dotplot.features) + RotatedAxis() & theme(axis.title.x = element_blank(),axis.title.y = element_blank())
		p
		
	#remaining cells have lower library sizes, but have more plasma cell and IgG expression, and will thus be annotated as IgG
		Idents(SO.integrated, cells = WhichCells(SO.integrated, idents = "plasmocyte"))     <- "IgG plasma cell" 
	
	#show plot of cells
		DimPlot(SO.integrated, label=T)
	
	#save annotations in annotations.list
		annotations.list <- list()
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents ="IgA plasma cell" )        ;names(annotations.list)[length(annotations.list)] <-   "IgA plasma cell"           
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents ="IgG plasma cell" )        ;names(annotations.list)[length(annotations.list)] <-   "IgG plasma cell"            
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents ="IgM plasma cell" )        ;names(annotations.list)[length(annotations.list)] <-   "IgM plasma cell"           
	
	#reset seurat object
		SO.integrated <- SO.backup
	
	#Inject new annotations in Seurat object
		for(i in seq_along(annotations.list)){
			Idents(SO.integrated, cells = annotations.list[[i]])  <- names(annotations.list)[[i]]
		}
	
	#save in meta.data$annotation
		SO.integrated[["annotation_2"]] <- Idents(SO.integrated)
		DimPlot(SO.integrated, label=T, label.size=5)

	#reorder idents and create DimPlots 
		Idents(SO.integrated) <- factor(Idents(SO.integrated), levels = c("IgA plasma cell","IgG plasma cell","IgM plasma cell","hepatocyte","hep/T/NK doublet","cholangiocyte","endothelial (1)","endothelial (2)","endothelial (3)","stromal","Tn/Tcm","T cytotox","Texh","T mito hi","MAIT/ILC (1)","MAIT/ILC (2)","NK CD56bright","NK CD56dim","B cell","mono CL","mono NCL","cDC1 (1)","cDC1 (2)","cDC2 (1)","cDC2 (2)","KC (1)","KC (2)","pDC","neutrophil","mast cell","cycling"))
		DimPlot(SO.integrated, label=T, label.size=5)
	
	# save 
		#saveRDS(SO.integrated, file="Narmada_pathseq_annotation_2.Rds")      #save
		#SO.integrated <- readRDS(file="/mnt/bctl/nvan0112/UZA/revision_JHEPR/data_files/scrnaseq_data/Narmada_pathseq_annotation_2.Rds")

	#remove cycling
		#SO.integrated <- subset(SO.integrated, idents="cycling", invert=TRUE)	

    #subset "no Sloss"
        table(SO.integrated@meta.data$SLoss)
        SO.integrated <- subset(SO.integrated, SLoss == "NoSLoss")

		#subset NAIVE viraemic patients (according to Table 1 Narmada et al)
		    table(SO.integrated@meta.data$orig.ident)
    	    SO.integrated <- subset(SO.integrated, orig.ident %in% c("CHB14","CHB9","CHB8","CHB7","CHB5","CHB4","CHB28","CHB24","CHB23","CHB3","CHB2"))
			table(SO.integrated@meta.data$orig.ident)

#plot
	#set colors
		n.clusters <- length(unique(Idents(SO.integrated)))
		if(n.clusters <= 10){
		  color.vector <- pal_npg("nrc")(n.clusters)
		}else{
		  #retrieve color codes
		  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
		    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
		    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		  }
		  color.vector <- c(pal_npg("nrc")(10),ggplotColours(n= (n.clusters-10)))
		}
		
	#UMAP representation plots
		p=DimPlot(SO.integrated, cols= color.vector, alpha=0.5, pt.size=0.5) + theme_void() + theme(legend.position = 'none')
		p
		if(save.plots==TRUE){  ggsave(p, file= "plots/DimPlot_HBV_viraemic_void_a5.jpg", width = 5, height=4)  }
		
		#labeled
			LabelClusters(DimPlot(SO.integrated, cols= color.vector.alpha2, pt.size=0.5) + theme_void() + theme(legend.position = 'none')   , id = "ident", color = "black", size = 5, repel = T,  box.padding = 1)

		
# Zoom in on coordinates of plasma cells
	SO.zoom <- subset(SO.integrated, idents=c("IgA plasma cell","IgM plasma cell","IgG plasma cell"))
		
		##cut everything beyond these UMAP coordinates
		#UMAP.coordinates.top    <- -0.5
		#UMAP.coordinates.bottom <- -6
		#UMAP.coordinates.left   <- -7
		#UMAP.coordinates.right  <- -3
		#
		#remove.cell.top    <- which(SO.zoom[["umap"]]@cell.embeddings[,2] > UMAP.coordinates.top) %>% names
		#remove.cell.bottom <- which(SO.zoom[["umap"]]@cell.embeddings[,2] < UMAP.coordinates.bottom) %>% names
		#remove.cell.left   <- which(SO.zoom[["umap"]]@cell.embeddings[,1] < UMAP.coordinates.left) %>% names
		#remove.cell.right  <- which(SO.zoom[["umap"]]@cell.embeddings[,1] > UMAP.coordinates.right) %>% names
		#
		#SO.zoom <- subset(SO.zoom, cells= unique(c(remove.cell.top,remove.cell.bottom,remove.cell.left,remove.cell.right)), invert=TRUE)
		
		DimPlot(SO.zoom, cols=color.vector[1:3])
		
		p=DimPlot(SO.zoom, cols=color.vector[1:3], pt.size=2, alpha=0.5) +theme_void() + NoLegend()
		p
		if(save.plots==TRUE){  ggsave(p, file= "plots/DimPlot_HBV_viraemic_plasmocyte_zoom.jpg", width = 3, height=1.5)  }
		
	# create figure dotplot B and plasma cell markers
		library(RColorBrewer)
		SO.subset <- subset(SO.integrated, idents=c("IgA plasma cell","IgM plasma cell","IgG plasma cell","B cell"))
		Idents(SO.subset) <- factor(Idents(SO.subset), levels = c("B cell","IgM plasma cell","IgG plasma cell","IgA plasma cell"))
		
		dotplot.features <- c("MZB1","TNFRSF17","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","IGHM","IGHD","BANK1","MS4A1","TNFRSF13C")
		
		p=DotPlot(SO.subset, cols = "RdYlBu",features = dotplot.features) + RotatedAxis() &  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
		p
		if(save.plots==TRUE){ggsave(p, file= "plots/DimPlot_HBV_viraemic_dotplot_sc_B_Plasma.jpg", width = 8, height=2.5)}
		if(save.plots==TRUE){ggsave(p, file= "plots/DimPlot_HBV_viraemic_dotplot_sc_B_Plasma_legend.jpg", width = 8, height=4)} #with legend



	# sessionInfo()
	# R version 4.4.1 (2024-06-14)
	# Platform: x86_64-pc-linux-gnu
	# Running under: Ubuntu 20.04.6 LTS

	# Matrix products: default
	# BLAS/LAPACK: FlexiBLAS OPENBLAS;  LAPACK version 3.11.0

	# locale:
	#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
	#  [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
	# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

	# time zone: Europe/Brussels
	# tzcode source: system (glibc)

	# attached base packages:
	# [1] stats     graphics  grDevices utils     datasets  methods   base     

	# other attached packages:
	# [1] ggsci_3.2.0        RColorBrewer_1.1-3 ggplot2_3.5.2      Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0           dplyr_1.1.4        httpgd_2.0.2      

	# loaded via a namespace (and not attached):
	#   [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.6            magrittr_2.0.3         RcppAnnoy_0.0.22      
	#   [7] matrixStats_1.5.0      ggridges_0.5.6         compiler_4.4.1         spatstat.geom_3.4-1    png_0.1-8              systemfonts_1.2.3     
	#  [13] vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1          crayon_1.5.3           pkgconfig_2.0.3        fastmap_1.2.0         
	#  [19] labeling_0.4.3         promises_1.3.3         ggbeeswarm_0.7.2       purrr_1.0.4            jsonlite_2.0.0         goftest_1.2-3         
	#  [25] later_1.4.2            spatstat.utils_3.1-4   irlba_2.3.5.1          parallel_4.4.1         cluster_2.1.6          R6_2.6.1              
	#  [31] ica_1.0-3              stringi_1.8.7          spatstat.data_3.1-6    reticulate_1.42.0      parallelly_1.45.0      spatstat.univar_3.1-3 
	#  [37] lmtest_0.9-40          scattermore_1.2        Rcpp_1.1.0             tensor_1.5.1           future.apply_1.20.0    zoo_1.8-14            
	#  [43] sctransform_0.4.2      httpuv_1.6.16          Matrix_1.7-0           splines_4.4.1          igraph_2.1.4           tidyselect_1.2.1      
	#  [49] dichromat_2.0-0.1      abind_1.4-8            spatstat.random_3.4-1  codetools_0.2-20       miniUI_0.1.2           spatstat.explore_3.4-3
	#  [55] listenv_0.9.1          lattice_0.22-6         tibble_3.3.0           plyr_1.8.9             withr_3.0.2            shiny_1.11.1          
	#  [61] ROCR_1.0-11            ggrastr_1.0.2          Rtsne_0.17             future_1.58.0          fastDummies_1.7.5      unigd_0.1.2           
	#  [67] survival_3.7-0         polyclip_1.10-7        fitdistrplus_1.2-3     pillar_1.10.2          KernSmooth_2.23-24     plotly_4.11.0         
	#  [73] generics_0.1.4         RcppHNSW_0.6.0         scales_1.4.0           globals_0.18.0         xtable_1.8-4           glue_1.8.0            
	#  [79] lazyeval_0.2.2         tools_4.4.1            data.table_1.17.6      RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2         
	#  [85] cowplot_1.1.3          grid_4.4.1             tidyr_1.3.1            nlme_3.1-165           patchwork_1.3.1        beeswarm_0.4.0        
	#  [91] vipor_0.4.7            cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-1            viridisLite_0.4.2      uwot_0.2.3            
	#  [97] gtable_0.3.6           digest_0.6.37          progressr_0.15.1       ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2          
	# [103] htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7             mime_0.13              MASS_7.3-61     				
		