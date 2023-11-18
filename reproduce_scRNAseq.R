#####################################################################################
#
#	Process liver scRNAseq data
#	Originally described by MacParland et al. 2018 Nature Communications
#	Raw data matrices provided by Brendan T. Innes
#
#	For reference, please cite MacParland et al. 2018 Nature Communications
#	https://doi.org/10.1038/s41467-018-06318-7
#
#####################################################################################

#note: the output of this analysis is available as an .rds file for R containing the seurat object on https://zenodo.org/records/10149895
#note: DOI 10.5281/zenodo.10149894

#set working directory
	setwd("C:/my_dir")
	
#load libraries
	library(Seurat)
	library(harmony)
	library(ggsci)
	library(ggplot2)
	
#1. load single cell data 
	#Create sample name vector and data directory vector
		sample.name.vector    <- c("macparland.donor1","macparland.donor2","macparland.donor3","macparland.donor4","macparland.donor5")
		data.dir.vector       <- c("P1TLH","P2TLH","P3TLH","P4TLH","P5TLH")
		SO.list <- list()
	
	#create seurat object list with selected samples, process them and add meta.data
		for(i in seq_along(sample.name.vector)){ 
			SO.list[[i]] <- Read10X(data.dir = eval(paste0("data_files/scRNAseq_data/",data.dir.vector[i])))
			SO.list[[i]] <- CreateSeuratObject(counts = SO.list[[i]], project = sample.name.vector[i], min.cells = 3, min.features = 200) 
			SO.list[[i]][["percent.mt"]]    <- PercentageFeatureSet(SO.list[[i]], pattern = "^MT-")
			SO.list[[i]] <- RenameCells(SO.list[[i]], add.cell.id = data.dir.vector[i])
			SO.list[[i]] <- subset(SO.list[[i]], subset = nFeature_RNA > 200 & percent.mt < 20 & nCount_RNA > 1500 & nCount_RNA < 99999)  
		}

	#integrate using Harmony harmony
		SO.integrated <- merge(SO.list[[1]], y = SO.list[-1])
		rm(SO.list)
	
		SO.integrated <- NormalizeData(SO.integrated)
		SO.integrated <- FindVariableFeatures(SO.integrated, nfeatures = 2000)
		SO.integrated <- ScaleData(SO.integrated)
		SO.integrated <- RunPCA(SO.integrated)
		SO.integrated <- RunHarmony(SO.integrated, group.by.vars = "orig.ident")
		SO.integrated <- RunUMAP(SO.integrated, reduction = "harmony", dims = 1:30)
		SO.integrated <- FindNeighbors(SO.integrated, reduction = "harmony", dims = 1:30) 
		SO.integrated <- FindClusters(SO.integrated, resolution = 2) 

#2. delete RBC (cluster 20) and reintegrate
	#delete red blood cell cluster
		SO.integrated <- subset(SO.integrated, idents = 20, invert = TRUE)

	#trim each sample to raw RNA counts only
		merge.vector <- as.character(unique(SO.integrated@meta.data$orig.ident))
		SO.merged.list <- list()
		for (i in seq_along(merge.vector)) {
			SO.merged.list[[i]] <- subset(SO.integrated, subset = orig.ident == merge.vector[i]) 
			DefaultAssay(SO.merged.list[[i]]) <- "RNA"
			SO.merged.list[[i]] <- DietSeurat(SO.merged.list[[i]] , assays = "RNA")
		}

	#reintegrate
		SO.integrated <- merge(SO.merged.list[[1]], y = SO.merged.list[-1])
		SO.integrated <- NormalizeData(SO.integrated)
		SO.integrated <- FindVariableFeatures(SO.integrated, nfeatures = 3000)
		SO.integrated <- ScaleData(SO.integrated)
		SO.integrated <- RunPCA(SO.integrated)
		SO.integrated <- RunHarmony(SO.integrated, group.by.vars = "orig.ident")
		SO.integrated <- RunUMAP(SO.integrated, reduction = "harmony", dims = 1:30)
		SO.integrated <- FindNeighbors(SO.integrated, reduction = "harmony", dims = 1:30) 
		SO.integrated <- FindClusters(SO.integrated, resolution = 2) 

#3. delete endothelial/kupffer cell doublets and reintegrate
	#delete endothelial/kupffer cell doublets
		SO.integrated <- subset(SO.integrated, idents = 25, invert = TRUE)

	#trim each sample to raw RNA counts only
		merge.vector <- as.character(unique(SO.integrated@meta.data$orig.ident))
		SO.merged.list <- list()
		for (i in seq_along(merge.vector)) {
			SO.merged.list[[i]] <- subset(SO.integrated, subset = orig.ident == merge.vector[i]) 
			DefaultAssay(SO.merged.list[[i]]) <- "RNA"
			SO.merged.list[[i]] <- DietSeurat(SO.merged.list[[i]] , assays = "RNA")
		}

	#reintegrate
		SO.integrated <- merge(SO.merged.list[[1]], y = SO.merged.list[-1])
		SO.integrated <- NormalizeData(SO.integrated)
		SO.integrated <- FindVariableFeatures(SO.integrated, nfeatures = 3000)
		SO.integrated <- ScaleData(SO.integrated)
		SO.integrated <- RunPCA(SO.integrated)
		SO.integrated <- RunHarmony(SO.integrated, group.by.vars = "orig.ident")
		SO.integrated <- RunUMAP(SO.integrated, reduction = "harmony", dims = 1:30)
		SO.integrated <- FindNeighbors(SO.integrated, reduction = "harmony", dims = 1:30) 
		SO.integrated <- FindClusters(SO.integrated, resolution = 2) 

SO.backup <- SO.integrated

#4. annotate
{
	annotations.list <- list()
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =0 )  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"     
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =1 )  ;names(annotations.list)[length(annotations.list)] <-   "monocyte"     
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =2 )  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"     
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =3 )  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"     
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =4 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"      
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =5 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"     
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =6 )  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"     
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =7 )  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =8 )  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56bright"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =9 )  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =10)  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =11)  ;names(annotations.list)[length(annotations.list)] <-   "cDC"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =12)  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =13)  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =14)  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56dim"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =15)  ;names(annotations.list)[length(annotations.list)] <-   "cycling"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =16)  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =17)  ;names(annotations.list)[length(annotations.list)] <-   "kupffer cell"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =18)  ;names(annotations.list)[length(annotations.list)] <-   "cycling"                 
	#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =19)  ;names(annotations.list)[length(annotations.list)] <-   "IgA/M plasma cell"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =20)  ;names(annotations.list)[length(annotations.list)] <-   "T cell"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =21)  ;names(annotations.list)[length(annotations.list)] <-   "IgG plasma cell"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =22)  ;names(annotations.list)[length(annotations.list)] <-   "cholangiocyte"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =23)  ;names(annotations.list)[length(annotations.list)] <-   "B cell"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =24)  ;names(annotations.list)[length(annotations.list)] <-   "cycling"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =25)  ;names(annotations.list)[length(annotations.list)] <-   "stromal cell"                 
	annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =26)  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"                 


	#divide plasma cell into IgA/M 
		SO.plasmaMA <- subset(SO.integrated, idents= "19")
			annotations.list[[length(annotations.list)+1]] <- WhichCells(subset(SO.plasmaMA, IGHM >= 4))   ;names(annotations.list)[length(annotations.list)] <-   "IgM plasma cell"     
			annotations.list[[length(annotations.list)+1]] <- WhichCells(subset(SO.plasmaMA, IGHM < 4))   ;names(annotations.list)[length(annotations.list)] <-   "IgA plasma cell"     
									
	# Inject new annotations
		for(i in seq_along(annotations.list)){
			Idents(SO.integrated, cells = annotations.list[[i]])  <- names(annotations.list)[[i]]
		}
		
	# save in meta.data$annotation
		SO.integrated[["annotation"]] <- Idents(SO.integrated)
}
	
#5. remove cycling
	SO.integrated <- subset(SO.integrated, idents = "cycling", invert = TRUE)

#6. reorder idents and create DimPlots 
	###del### Idents(SO.integrated) <- factor(Idents(SO.integrated), levels = c("hepatocyte","cholangiocyte","endothelial cell","stromal cell","T cell","NK CD56bright","NK CD56dim","monocyte","cDC","kupffer cell","B cell","IgM plasma cell","IgG plasma cell","IgA plasma cell"))
	###del### Idents(SO.integrated) <- factor(Idents(SO.integrated), levels = c("B cell","IgA plasma cell","IgG plasma cell","IgM plasma cell","hepatocyte","cholangiocyte","endothelial cell","stromal cell","T cell","NK CD56bright","NK CD56dim","monocyte","cDC","kupffer cell"))
	Idents(SO.integrated) <- factor(Idents(SO.integrated), levels = c("IgA plasma cell","IgG plasma cell","IgM plasma cell","B cell","cholangiocyte","stromal cell","hepatocyte","endothelial cell","T cell","NK CD56bright","NK CD56dim","monocyte","cDC","kupffer cell"))

#7. save 
	#saveRDS(SO.integrated, file="data_files/scRNAseq_data/liver_scRNAseq_IgAGM.Rds")      #save

#8. plot
	#rework colors
	color.vector <- pal_npg("nrc")(10)
	color.vector <- c(color.vector,"brown","khaki","orange3","orchid2")
	color.vector.alpha2 <- alpha(color.vector,0.2)	
	color.vector.alpha5 <- alpha(color.vector,0.5)	
 
	#UMAP representation plots
		DimPlot(SO.integrated, cols= color.vector, pt.size=1)
		DimPlot(SO.integrated, cols= color.vector, pt.size=0.5) + theme_void() + theme(legend.position = 'none')
		DimPlot(SO.integrated, cols= color.vector.alpha2, pt.size=0.5) + theme_void() + theme(legend.position = 'none')
		

	#DotPlots
		dotplot.markers <- c("MZB1","TNFRSF17","MS4A1","SOX9","KRT19","MYL9","COL1A1","TTR","ALB","CD34","PECAM1","CD2","CD3D","NCAM1","XCL1","XCL2","CD160","GNLY","GZMB","FCGR3A","FCN1","VCAN","HLA-DQA1","HLA-DPB1","HLA-DMA","HLA-DRB1","CLEC10A","MARCO","CD5L")
		DotPlot(SO.integrated, cols = "RdYlBu",features = dotplot.markers) + RotatedAxis() & theme(axis.title.x = element_blank(),axis.title.y = element_blank())
		
	#Cluster proportion table
		clustertable <- table(Idents(SO.integrated), SO.integrated@meta.data$orig.ident) #color by dataset (orig.ident)
		plot.orig.ident <- ggplot(as.data.frame(clustertable), aes(fill=Var2, y=Freq, x=Var1)) + 
			geom_bar(position="fill", stat="identity")
		
		denominator <-  as.numeric(colSums(clustertable))	
		clustertable.normalized <- clustertable
		for(i in 1:ncol(clustertable)){
			clustertable.normalized[,i] <- clustertable[,i]/denominator[i]
		}
		
		p=ggplot(as.data.frame(clustertable.normalized), aes(fill=Var2, y=Freq, x=Var1)) + 
			geom_bar(position="fill", stat="identity") + 
			theme_classic() +
			xlab("") +
			ylab("relative contribution") +
			theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.25, size = 20, color="black")) +
			theme(axis.text.y = element_text(size = 20, color="black")) +
			theme(axis.title.y = element_text(size = 20)) +
			scale_y_continuous(expand = c(0, 0), limits=c(0,1.05), breaks=c(0,0.25,0.5,0.75,1)) +
			theme(legend.text=element_text(size=20)) +
			theme(legend.title=element_text(size=20)) +
			guides(fill=guide_legend("donor sample"))
		p		
		
		
	# Zoom in on coordinates of plasma cells
		library(dplyr)
		#cut everything beyond these UMAP coordinates
			UMAP.coordinates.top    <- 2
			UMAP.coordinates.bottom <- -5
			UMAP.coordinates.left   <- -10
			UMAP.coordinates.right  <- -6.75
		
			remove.cell.top    <- which(SO.integrated[["umap"]]@cell.embeddings[,2] > UMAP.coordinates.top) %>% names
			remove.cell.bottom <- which(SO.integrated[["umap"]]@cell.embeddings[,2] < UMAP.coordinates.bottom) %>% names
			remove.cell.left   <- which(SO.integrated[["umap"]]@cell.embeddings[,1] < UMAP.coordinates.left) %>% names
			remove.cell.right  <- which(SO.integrated[["umap"]]@cell.embeddings[,1] > UMAP.coordinates.right) %>% names
			
			SO.zoom <- subset(SO.integrated, cells= unique(c(remove.cell.top,remove.cell.bottom,remove.cell.left,remove.cell.right)), invert=TRUE)
			DimPlot(SO.zoom, cols=color.vector[1:3])
		
		DimPlot(SO.zoom, cols=color.vector.alpha5[1:3], pt.size=3) +theme_void() + NoLegend()

	# create figure dotplot B and plasma cell markers
	library(RColorBrewer)
		SO.subset <- subset(SO.integrated, idents=c("IgA plasma cell","IgM plasma cell","IgG plasma cell","B cell"))
		Idents(SO.subset) <- factor(Idents(SO.subset), levels = c("B cell","IgG plasma cell","IgM plasma cell","IgA plasma cell"))
	
		dotplot.features <- c("MZB1","DERL3","TNFRSF17","POU2AF1","IGHA1","IGHA2","IGHM","IGHG1","IGHG2","IGHG3","IGHG4","BANK1","MS4A1","TNFRSF13C","CD22")
		
		DotPlot(SO.subset, cols = "RdYlBu",features = dotplot.features) + RotatedAxis() &  theme(axis.title.x = element_blank(),axis.title.y = element_blank())

	
##############################
#  SESSIONINFO()
##############################
#sessionInfo()
#
#	R version 4.0.3 (2020-10-10)
#	Platform: x86_64-w64-mingw32/x64 (64-bit)
#	Running under: Windows 10 x64 (build 19044)
#	
#	Matrix products: default
#	
#	locale:
#	[1] LC_COLLATE=Dutch_Belgium.1252  LC_CTYPE=Dutch_Belgium.1252    LC_MONETARY=Dutch_Belgium.1252 LC_NUMERIC=C                   LC_TIME=Dutch_Belgium.1252    
#	
#	attached base packages:
#	[1] stats     graphics  grDevices utils     datasets  methods   base     
#	
#	other attached packages:
#	[1] ggplot2_3.4.0      ggsci_2.9          harmony_0.1.0      Rcpp_1.0.7         SeuratObject_4.0.4 Seurat_4.0.2      
#	
#	loaded via a namespace (and not attached):
#	  [1] nlme_3.1-149          matrixStats_0.58.0    spatstat.sparse_2.0-0 RcppAnnoy_0.0.18      RColorBrewer_1.1-2    httr_1.4.2            sctransform_0.3.2     tools_4.0.3          
#	  [9] utf8_1.2.1            R6_2.5.0              irlba_2.3.3           rpart_4.1-15          KernSmooth_2.23-17    uwot_0.1.10           mgcv_1.8-33           DBI_1.1.1            
#	 [17] lazyeval_0.2.2        colorspace_2.0-0      withr_2.5.0           tidyselect_1.1.0      gridExtra_2.3         compiler_4.0.3        cli_3.4.1             plotly_4.9.3         
#	 [25] labeling_0.4.2        scales_1.2.1          lmtest_0.9-38         spatstat.data_2.1-0   ggridges_0.5.3        pbapply_1.4-3         goftest_1.2-2         stringr_1.4.0        
#	 [33] digest_0.6.27         spatstat.utils_2.3-1  pkgconfig_2.0.3       htmltools_0.5.2       parallelly_1.32.0     fastmap_1.1.0         htmlwidgets_1.5.3     rlang_1.0.6          
#	 [41] shiny_1.6.0           farver_2.1.0          generics_0.1.0        zoo_1.8-9             jsonlite_1.7.2        ica_1.0-2             dplyr_1.0.5           magrittr_2.0.1       
#	 [49] patchwork_1.1.1       Matrix_1.3-4          munsell_0.5.0         fansi_0.4.2           abind_1.4-5           reticulate_1.18       lifecycle_1.0.3       stringi_1.5.3        
#	 [57] MASS_7.3-53           Rtsne_0.15            plyr_1.8.6            grid_4.0.3            parallel_4.0.3        listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1        
#	 [65] crayon_1.4.1          miniUI_0.1.1.1        deldir_1.0-6          lattice_0.20-41       cowplot_1.1.1         splines_4.0.3         tensor_1.5            pillar_1.6.0         
#	 [73] igraph_1.2.6          spatstat.geom_2.4-0   future.apply_1.7.0    reshape2_1.4.4        codetools_0.2-16      leiden_0.3.7          glue_1.4.2            data.table_1.14.0    
#	 [81] png_0.1-7             vctrs_0.5.1           httpuv_1.5.5          gtable_0.3.0          RANN_2.6.1            purrr_0.3.4           spatstat.core_2.0-0   polyclip_1.10-0      
#	 [89] tidyr_1.1.3           scattermore_0.7       future_1.26.1         assertthat_0.2.1      mime_0.10             xtable_1.8-4          RSpectra_0.16-0       later_1.1.0.1        
#	 [97] survival_3.2-7        viridisLite_0.3.0     tibble_3.1.0          cluster_2.1.0         globals_0.15.1        fitdistrplus_1.1-3    ellipsis_0.3.2        ROCR_1.0-11          


