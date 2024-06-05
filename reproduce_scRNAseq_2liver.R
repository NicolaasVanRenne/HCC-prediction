#####################################################################################
#
#	Process liver scRNAseq data of two liver scRNAseq data sets
#	Originally described by Ramachandran et al. 2019 Nature 
#	Originally described by Nkongolo et al. 2023 J Clin Invest
#
#	For reference, please cite 
#	Ramachandran: DOI: 10.1038/s41586-019-1631-3
#	Nkongolo: DOI: 10.1172/JCI158903
#
#####################################################################################
#note: the output of this analysis is available as an .rds file for R containing the seurat object on https://zenodo.org/records/11487048
#note: DOI 10.5281/zenodo.11487047

#set working directory
	setwd("C:/my_dir")
	
#load packages
	library(Seurat)
	library(harmony)
	library(ggsci)
	library(ggplot2)
	library(RColorBrewer)
	library(patchwork)
	library(dplyr)
	library(ggsci)

#set memory
#options(future.globals.maxSize = 32000 * 1024^2)

###del### #set seed
###del### 	seed = 18101983
###del### 	set.seed(seed)

#load scRNAseq data
	#load settings
		feat.min.cells = 3 #minimal cell expressing threshold. Genes expressed in fewer cells will be omitted
		min.nFeature_RNA = 200 #minimal genes per cell
		min.nCount_RNA = 1500 #RNA count lower library
		max.nCount_RNA = 99999 #RNA count upper limit
		set.mito = 10 #mitochondrial content limit 
{
	SO.list <- list()
	SO.list.1 <- list()
	SO.list.2 <- list()
}

{	
	#########################################################################################################################################################
	# Load the Ramachandran datasets,
	# Ramachandran P, Nature 2019, resolving the fibrotic niche of human liver cirrhosis at single-cell level, GEO Series GSE136103
	#########################################################################################################################################################
	
	### data download instructions ###  The data (human liver only; ignore mouse and blood) can be downloaded from NCBI GEO GSE136103, store it in a folder called 'data_directory'
	### data download instructions ###  Each sample should be stored in a separate folder containing the folder structure in the data.dir.vector below
	### data download instructions ###  e.g; human.healthy1.CD45pos should be stored in subfolder 'Human_healthy1_CD45' and contains the following files:
	### data download instructions ###  Supplementary file	Size	Download	File type/resource
	### data download instructions ###  GSM4041150_healthy1_cd45+_barcodes.tsv.gz	8.4 Kb	(ftp)(http)	TSV
	### data download instructions ###  GSM4041150_healthy1_cd45+_genes.tsv.gz	258.6 Kb	(ftp)(http)	TSV
	### data download instructions ###  GSM4041150_healthy1_cd45+_matrix.mtx.gz	6.3 Mb	(ftp)(http)	MTX
	
	
	#Create sample name vector and data directory vector
		sample.name.vector    <- c("human.healthy1.CD45pos","human.healthy2.CD45pos","human.healthy3.CD45pos","human.healthy4.CD45pos","human.healthy5.CD45pos","human.healthy1.CD45negA","human.healthy1.CD45negB","human.healthy2.CD45neg","human.healthy3.CD45negA","human.healthy3.CD45negB","human.healthy4.CD45neg","human.cirrhotic1.CD45pos","human.cirrhotic2.CD45pos","human.cirrhotic3.CD45pos","human.cirrhotic4.CD45pos","human.cirrhotic5.CD45pos","human.cirrhotic1.CD45negA","human.cirrhotic1.CD45negB","human.cirrhotic2.CD45neg","human.cirrhotic3.CD45neg")
		data.dir.vector       <- c("Human_healthy1_CD45+","Human_healthy2_CD45+","Human_healthy3_CD45+","Human_healthy4_CD45+","Human_healthy5_CD45+","Human_healthy1_CD45-A","Human_healthy1_CD45-B","Human_healthy2_CD45-","Human_healthy3_CD45-A","Human_healthy3_CD45-B","Human_healthy4_CD45-","Human_cirrhotic1_CD45+","Human_cirrhotic2_CD45+","Human_cirrhotic3_CD45+","Human_cirrhotic4_CD45+","Human_cirrhotic5_CD45+","Human_cirrhotic1_CD45-A","Human_cirrhotic1_CD45-B","Human_cirrhotic2_CD45-","Human_cirrhotic3_CD45-")
	
		sample.name.vector.ramachandran <- sample.name.vector #save for later
		
	#create seurat object list with selected samples, process them and add meta.data
		for(i in seq_along(sample.name.vector)){ 
			SO.list.1[[i]] <- Read10X(data.dir = eval(paste0("data_directory/",data.dir.vector[i])))
			SO.list.1[[i]] <- CreateSeuratObject(counts = SO.list.1[[i]], project = sample.name.vector[i], min.cells = feat.min.cells, min.features = min.nFeature_RNA) 
			if(is.element(i,1:24))  {SO.list.1[[i]][["percent.mt"]] <- PercentageFeatureSet(SO.list.1[[i]], pattern = "^MT-")}
			if(is.element(i,c(1:11)))   {SO.list.1[[i]][["status"]] <- "healthy"}
			if(is.element(i,c(12:24)))  {SO.list.1[[i]][["status"]] <- "cirrhotic"}
			if(is.element(i,c(1,6,7)))      {SO.list.1[[i]][["patient"]] <- "Ramach_P01"}
			if(is.element(i,c(2,8)))        {SO.list.1[[i]][["patient"]] <- "Ramach_P02"}
			if(is.element(i,c(3,9,10)))     {SO.list.1[[i]][["patient"]] <- "Ramach_P03"}
			if(is.element(i,c(4,11)))       {SO.list.1[[i]][["patient"]] <- "Ramach_P04"}
			if(is.element(i,c(5)))          {SO.list.1[[i]][["patient"]] <- "Ramach_P05"}
			if(is.element(i,c(12,17,18)))   {SO.list.1[[i]][["patient"]] <- "Ramach_P06"}
			if(is.element(i,c(13,19)))      {SO.list.1[[i]][["patient"]] <- "Ramach_P07"}
			if(is.element(i,c(14,20)))      {SO.list.1[[i]][["patient"]] <- "Ramach_P08"}
			if(is.element(i,15))            {SO.list.1[[i]][["patient"]] <- "Ramach_P09"}
			if(is.element(i,16))            {SO.list.1[[i]][["patient"]] <- "Ramach_P10"}
			SO.list.1[[i]] <- RenameCells(SO.list.1[[i]], add.cell.id = data.dir.vector[i])
			SO.list.1[[i]] <- subset(SO.list.1[[i]], subset = nFeature_RNA > min.nFeature_RNA & percent.mt < set.mito & nCount_RNA > min.nCount_RNA & nCount_RNA < max.nCount_RNA)  
		}
	}
	
	{
	#########################################################################################################################################################
	# Load the nkongolo2023 datasets,
	# nkongolo S, J Clin Invest 2023, Longitudinal liver sampling in patients with chronic hepatitis B starting antiviral therapy reveals hepatotoxic CD8+ T cells
	#########################################################################################################################################################
	
	### data download instructions ### 	# The data can be downloaded from NCBI GEO GSE216314, store it in a folder called 'data_directory'
	### data download instructions ###   # Supplementary file	Size	Download	File type/resource
	### data download instructions ###   # GSE216314_Metadata.txt.gz	234.5 Kb	(ftp)(http)	TXT
	### data download instructions ###   # GSE216314_barcodes.tsv.gz	222.3 Kb	(ftp)(http)	TSV
	### data download instructions ###   # GSE216314_features.tsv.gz	310.8 Kb	(ftp)(http)	TSV
	### data download instructions ###   # GSE216314_matrix.mtx.gz	222.6 Mb	(ftp)(http)	MTX
	
	### run only first time!!!### 	#First we need to separate the samples for integration, because everything is in one big matrix
	### run only first time!!!### 	#load data
	### run only first time!!!### 		nkongolo.data <- Read10X(data.dir = "data_directory", gene.column=1)
	### run only first time!!!### 	
	### run only first time!!!### 	#identify donor samples and split dataset
	### run only first time!!!### 		barcodes <- colnames(nkongolo.data)
	### run only first time!!!### 	
	### run only first time!!!### 		metadata <- read.delim(file="data_directory/GSE216314_Metadata.txt.gz", header=T)
	### run only first time!!!### 	
	### run only first time!!!### 		samples <- substr(barcodes, 1,8)
	### run only first time!!!### 		samples <- gsub("WK0_","WK00",samples)
	### run only first time!!!### 	
	### run only first time!!!### 		metadata$samples <- samples
	### run only first time!!!### 	
	### run only first time!!!### 	
	### run only first time!!!### 	P004WK00 <- metadata$CellID[which(metadata$samples == "004_WK00")]
	### run only first time!!!### 	P004WK12 <- metadata$CellID[which(metadata$samples == "004_WK12")]
	### run only first time!!!### 	P004WK24 <- metadata$CellID[which(metadata$samples == "004_WK24")]
	### run only first time!!!### 	P008WK00 <- metadata$CellID[which(metadata$samples == "008_WK00")]
	### run only first time!!!### 	P008WK12 <- metadata$CellID[which(metadata$samples == "008_WK12")]
	### run only first time!!!### 	P008WK24 <- metadata$CellID[which(metadata$samples == "008_WK24")]
	### run only first time!!!### 	P009WK00 <- metadata$CellID[which(metadata$samples == "009_WK00")]
	### run only first time!!!### 	P009WK12 <- metadata$CellID[which(metadata$samples == "009_WK12")]
	### run only first time!!!### 	P009WK24 <- metadata$CellID[which(metadata$samples == "009_WK24")]
	### run only first time!!!### 	P010WK00 <- metadata$CellID[which(metadata$samples == "010_WK00")]
	### run only first time!!!### 	P010WK12 <- metadata$CellID[which(metadata$samples == "010_WK12")]
	### run only first time!!!### 	P010WK24 <- metadata$CellID[which(metadata$samples == "010_WK24")]
	### run only first time!!!### 	P011WK00 <- metadata$CellID[which(metadata$samples == "011_WK00")]
	### run only first time!!!### 	P011WK12 <- metadata$CellID[which(metadata$samples == "011_WK12")]
	### run only first time!!!### 	P011WK24 <- metadata$CellID[which(metadata$samples == "011_WK24")]
	### run only first time!!!### 	
	### run only first time!!!### 	P004WK00.data <- nkongolo.data[,P004WK00]   ;  saveRDS(P004WK00.data , file="data_directory/nkongolo2023_P004WK00.rds")
	### run only first time!!!### 	P004WK12.data <- nkongolo.data[,P004WK12]   ;  saveRDS(P004WK12.data , file="data_directory/nkongolo2023_P004WK12.rds")
	### run only first time!!!### 	P004WK24.data <- nkongolo.data[,P004WK24]   ;  saveRDS(P004WK24.data , file="data_directory/nkongolo2023_P004WK24.rds")
	### run only first time!!!### 	P008WK00.data <- nkongolo.data[,P008WK00]   ;  saveRDS(P008WK00.data , file="data_directory/nkongolo2023_P008WK00.rds")
	### run only first time!!!### 	P008WK12.data <- nkongolo.data[,P008WK12]   ;  saveRDS(P008WK12.data , file="data_directory/nkongolo2023_P008WK12.rds")
	### run only first time!!!### 	P008WK24.data <- nkongolo.data[,P008WK24]   ;  saveRDS(P008WK24.data , file="data_directory/nkongolo2023_P008WK24.rds")
	### run only first time!!!### 	P009WK00.data <- nkongolo.data[,P009WK00]   ;  saveRDS(P009WK00.data , file="data_directory/nkongolo2023_P009WK00.rds")
	### run only first time!!!### 	P009WK12.data <- nkongolo.data[,P009WK12]   ;  saveRDS(P009WK12.data , file="data_directory/nkongolo2023_P009WK12.rds")
	### run only first time!!!### 	P009WK24.data <- nkongolo.data[,P009WK24]   ;  saveRDS(P009WK24.data , file="data_directory/nkongolo2023_P009WK24.rds")
	### run only first time!!!### 	P010WK00.data <- nkongolo.data[,P010WK00]   ;  saveRDS(P010WK00.data , file="data_directory/nkongolo2023_P010WK00.rds")
	### run only first time!!!### 	P010WK12.data <- nkongolo.data[,P010WK12]   ;  saveRDS(P010WK12.data , file="data_directory/nkongolo2023_P010WK12.rds")
	### run only first time!!!### 	P010WK24.data <- nkongolo.data[,P010WK24]   ;  saveRDS(P010WK24.data , file="data_directory/nkongolo2023_P010WK24.rds")
	### run only first time!!!### 	P011WK00.data <- nkongolo.data[,P011WK00]   ;  saveRDS(P011WK00.data , file="data_directory/nkongolo2023_P011WK00.rds")
	### run only first time!!!### 	P011WK12.data <- nkongolo.data[,P011WK12]   ;  saveRDS(P011WK12.data , file="data_directory/nkongolo2023_P011WK12.rds")
	### run only first time!!!### 	P011WK24.data <- nkongolo.data[,P011WK24]   ;  saveRDS(P011WK24.data , file="data_directory/nkongolo2023_P011WK24.rds")
	
	
	#Create sample name vector and data directory vector
		sample.name.vector    <- c("nkongolo2023_P004WK00","nkongolo2023_P004WK12","nkongolo2023_P004WK24","nkongolo2023_P008WK00","nkongolo2023_P008WK12","nkongolo2023_P008WK24","nkongolo2023_P009WK00","nkongolo2023_P009WK12","nkongolo2023_P009WK24","nkongolo2023_P010WK00","nkongolo2023_P010WK12","nkongolo2023_P010WK24","nkongolo2023_P011WK00","nkongolo2023_P011WK12","nkongolo2023_P011WK24")
		sample.name.vector.nkongolo <- sample.name.vector #save for later
		
	#create seurat object list with selected samples, process them and add meta.data
		for(i in seq_along(sample.name.vector)){ 
			SO.list.2[[i]] <- readRDS(file = eval(paste0("data_directory/",sample.name.vector[i],".rds")))
			SO.list.2[[i]] <- CreateSeuratObject(counts = SO.list.2[[i]], project = sample.name.vector[i], min.cells = feat.min.cells, min.features = min.nFeature_RNA) 
			SO.list.2[[i]][["percent.mt"]]    <- PercentageFeatureSet(SO.list.2[[i]], pattern = "^MT-")
			SO.list.2[[i]][["tissue"]] <- "liver" 
			SO.list.2[[i]] <- subset(SO.list.2[[i]], subset = nFeature_RNA > min.nFeature_RNA & percent.mt < set.mito & nCount_RNA > min.nCount_RNA & nCount_RNA < max.nCount_RNA)  
	
			#set orig.ident
			if(is.element(i,1 ))     {SO.list.2[[i]][["orig.ident"]] <- "P004WK00"}
			if(is.element(i,2 ))     {SO.list.2[[i]][["orig.ident"]] <- "P004WK12"}
			if(is.element(i,3 ))     {SO.list.2[[i]][["orig.ident"]] <- "P004WK24"}
			if(is.element(i,4 ))     {SO.list.2[[i]][["orig.ident"]] <- "P008WK00"}
			if(is.element(i,5 ))     {SO.list.2[[i]][["orig.ident"]] <- "P008WK12"}
			if(is.element(i,6 ))     {SO.list.2[[i]][["orig.ident"]] <- "P008WK24"}
			if(is.element(i,7 ))     {SO.list.2[[i]][["orig.ident"]] <- "P009WK00"}
			if(is.element(i,8 ))     {SO.list.2[[i]][["orig.ident"]] <- "P009WK12"}
			if(is.element(i,9 ))     {SO.list.2[[i]][["orig.ident"]] <- "P009WK24"}
			if(is.element(i,10))     {SO.list.2[[i]][["orig.ident"]] <- "P010WK00"}
			if(is.element(i,11))     {SO.list.2[[i]][["orig.ident"]] <- "P010WK12"}
			if(is.element(i,12))     {SO.list.2[[i]][["orig.ident"]] <- "P010WK24"}
			if(is.element(i,13))     {SO.list.2[[i]][["orig.ident"]] <- "P011WK00"}
			if(is.element(i,14))     {SO.list.2[[i]][["orig.ident"]] <- "P011WK12"}
			if(is.element(i,15))     {SO.list.2[[i]][["orig.ident"]] <- "P011WK24"}
			
			#add metadata: patient
			if(is.element(i,1 ))     {SO.list.2[[i]][["patient"]] <- "P004"}
			if(is.element(i,2 ))     {SO.list.2[[i]][["patient"]] <- "P004"}
			if(is.element(i,3 ))     {SO.list.2[[i]][["patient"]] <- "P004"}
			if(is.element(i,4 ))     {SO.list.2[[i]][["patient"]] <- "P008"}
			if(is.element(i,5 ))     {SO.list.2[[i]][["patient"]] <- "P008"}
			if(is.element(i,6 ))     {SO.list.2[[i]][["patient"]] <- "P008"}
			if(is.element(i,7 ))     {SO.list.2[[i]][["patient"]] <- "P009"}
			if(is.element(i,8 ))     {SO.list.2[[i]][["patient"]] <- "P009"}
			if(is.element(i,9 ))     {SO.list.2[[i]][["patient"]] <- "P009"}
			if(is.element(i,10))     {SO.list.2[[i]][["patient"]] <- "P010"}
			if(is.element(i,11))     {SO.list.2[[i]][["patient"]] <- "P010"}
			if(is.element(i,12))     {SO.list.2[[i]][["patient"]] <- "P010"}
			if(is.element(i,13))     {SO.list.2[[i]][["patient"]] <- "P011"}
			if(is.element(i,14))     {SO.list.2[[i]][["patient"]] <- "P011"}
			if(is.element(i,15))     {SO.list.2[[i]][["patient"]] <- "P011"}
	
			#add metadata: week
			if(is.element(i,1 ))     {SO.list.2[[i]][["week"]] <- "week0"}
			if(is.element(i,2 ))     {SO.list.2[[i]][["week"]] <- "week12"}
			if(is.element(i,3 ))     {SO.list.2[[i]][["week"]] <- "week24"}
			if(is.element(i,4 ))     {SO.list.2[[i]][["week"]] <- "week0"}
			if(is.element(i,5 ))     {SO.list.2[[i]][["week"]] <- "week12"}
			if(is.element(i,6 ))     {SO.list.2[[i]][["week"]] <- "week24"}
			if(is.element(i,7 ))     {SO.list.2[[i]][["week"]] <- "week0"}
			if(is.element(i,8 ))     {SO.list.2[[i]][["week"]] <- "week12"}
			if(is.element(i,9 ))     {SO.list.2[[i]][["week"]] <- "week24"}
			if(is.element(i,10))     {SO.list.2[[i]][["week"]] <- "week0"}
			if(is.element(i,11))     {SO.list.2[[i]][["week"]] <- "week12"}
			if(is.element(i,12))     {SO.list.2[[i]][["week"]] <- "week24"}
			if(is.element(i,13))     {SO.list.2[[i]][["week"]] <- "week0"}
			if(is.element(i,14))     {SO.list.2[[i]][["week"]] <- "week12"}
			if(is.element(i,15))     {SO.list.2[[i]][["week"]] <- "week24"}
		}
	}
	
		#merge lists of seurat objects and integrate using harmony
			SO.list <- c(SO.list.2,SO.list.1)
			rm(SO.list.1, SO.list.2);gc()
	
	
			SO.integrated <- merge(SO.list[[1]], y = SO.list[-1])
			rm(SO.list);gc()
		
			SO.integrated <- NormalizeData(SO.integrated)
			SO.integrated <- FindVariableFeatures(SO.integrated, nfeatures = 2000)
			SO.integrated <- ScaleData(SO.integrated)
			SO.integrated <- RunPCA(SO.integrated)
			SO.integrated <- RunHarmony(SO.integrated, group.by.vars = "orig.ident")
			SO.integrated <- RunUMAP(SO.integrated, reduction = "harmony", dims = 1:30)
			SO.integrated <- FindNeighbors(SO.integrated, reduction = "harmony", dims = 1:30) 
			SO.integrated <- FindClusters(SO.integrated, resolution = 2) 
		
	
			DimPlot(SO.integrated, label=T, label.size=7)

	#1. set annotations 
		annotations.list <- list()
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =0 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =1 )  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56bright"     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =2 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =3 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =4 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =5 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =6 )  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56dim"     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =7 )  ;names(annotations.list)[length(annotations.list)] <-   "myeloid cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =8 )  ;names(annotations.list)[length(annotations.list)] <-   "myeloid cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =9 )  ;names(annotations.list)[length(annotations.list)] <-   "T cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =10)  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56dim"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =11)  ;names(annotations.list)[length(annotations.list)] <-   "B cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =12)  ;names(annotations.list)[length(annotations.list)] <-   "myeloid cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =13)  ;names(annotations.list)[length(annotations.list)] <-   "T cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =14)  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56bright"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =15)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =16)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =17)  ;names(annotations.list)[length(annotations.list)] <-   "cholangiocyte"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =18)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =19)  ;names(annotations.list)[length(annotations.list)] <-   "myeloid cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =20)  ;names(annotations.list)[length(annotations.list)] <-   "stromal cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =21)  ;names(annotations.list)[length(annotations.list)] <-   "myeloid cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =22)  ;names(annotations.list)[length(annotations.list)] <-   "cycling"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =23)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =24)  ;names(annotations.list)[length(annotations.list)] <-   "T cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =25)  ;names(annotations.list)[length(annotations.list)] <-   "25"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =26)  ;names(annotations.list)[length(annotations.list)] <-   "T cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =27)  ;names(annotations.list)[length(annotations.list)] <-   "myeloid cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =28)  ;names(annotations.list)[length(annotations.list)] <-   "myeloid cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =29)  ;names(annotations.list)[length(annotations.list)] <-   "stromal cell"                 
		#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =30)  ;names(annotations.list)[length(annotations.list)] <-   "epithelial cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =31)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =32)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =33)  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56dim"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =34)  ;names(annotations.list)[length(annotations.list)] <-   "cycling"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =35)  ;names(annotations.list)[length(annotations.list)] <-   "T cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =36)  ;names(annotations.list)[length(annotations.list)] <-   "T cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =37)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =38)  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56bright"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =39)  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56bright"                 
                                                       
	#2a. Subset and annotate plasma cells
	
		SO.integrated -> SO.backup
		
		SO.integrated <- subset(SO.integrated, idents=25)
		
		#rename plasma cells according to heavy chain expression
			Idents(SO.integrated, cells = WhichCells(SO.integrated, expression = `IGHG1` > 4)) <- "IgG plasma cell"
			Idents(SO.integrated, cells = WhichCells(SO.integrated,  expression = `IGHG2` > 4)) <- "IgG plasma cell"
			Idents(SO.integrated, cells = WhichCells(SO.integrated,  expression = `IGHG3` > 4)) <- "IgG plasma cell"
			Idents(SO.integrated, cells = WhichCells(SO.integrated,  expression = `IGHG4` > 4)) <- "IgG plasma cell"
			Idents(SO.integrated, cells = WhichCells(SO.integrated,  expression = `IGHA1` > 4)) <- "IgA plasma cell"
			Idents(SO.integrated, cells = WhichCells(SO.integrated,  expression = `IGHA2` > 4)) <- "IgA plasma cell"
			Idents(SO.integrated, cells = WhichCells(SO.integrated,  expression = `IGHM` > 4))  <- "IgM plasma cell"
			Idents(SO.integrated, cells = WhichCells(SO.integrated, idents = 25)) <- "IgD plasma cell"
			
			DimPlot(SO.integrated, label=T)
		
		#save annotations in annotations.list
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents ="IgA plasma cell" )        ;names(annotations.list)[length(annotations.list)] <-   "IgA plasma cell"           
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents ="IgG plasma cell" )        ;names(annotations.list)[length(annotations.list)] <-   "IgG plasma cell"            
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents ="IgM plasma cell" )        ;names(annotations.list)[length(annotations.list)] <-   "IgM plasma cell"           
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents ="IgD plasma cell" )        ;names(annotations.list)[length(annotations.list)] <-   "IgD plasma cell"            
	
	
		#reset seurat object
			SO.integrated <- SO.backup

	#2b. Subset and annotate epithelial cells
		#SO.integrated <- subset(SO.integrated, idents=c(17,30))
		SO.integrated <- subset(SO.integrated, idents=c(30))
	
		#change granularity
		{
			#1. set cluster resolution
			FindCluster.resolution = 0.4  
				
			#2. run algorithm
			DefaultAssay(SO.integrated) <- "RNA"
			
			if(length(SO.integrated@graphs) == 0){ #if graph is not present
						SO.integrated <- FindNeighbors(SO.integrated, reduction = "pca", dims = 1:30)
					}
			SO.integrated <- FindClusters(SO.integrated, resolution = FindCluster.resolution) #hoe hoger res, hoe meer clusters
			
			# plot
			plot.orig.ident <- DimPlot(SO.integrated, group.by = "orig.ident") + ggtitle(paste("change granularity"))
			plot.ident   <- DimPlot(SO.integrated,  label=TRUE, label.size=7)   + ggtitle(paste0("FindCluster.resolution = ",FindCluster.resolution))   
			wrap_plots(list(plot.orig.ident, plot.ident))
		}


		#add annotations to list
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =0 )  ;names(annotations.list)[length(annotations.list)] <-   "cholangiocyte"
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =1 )  ;names(annotations.list)[length(annotations.list)] <-   "cholangiocyte"
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =2 )  ;names(annotations.list)[length(annotations.list)] <-   "cholangiocyte"
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =3 )  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =4 )  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =5 )  ;names(annotations.list)[length(annotations.list)] <-   "cholangiocyte"
		

		#reset seurat object
		SO.integrated <- SO.backup	

	#3. Inject new annotations in Seurat object
		for(i in seq_along(annotations.list)){
			Idents(SO.integrated, cells = annotations.list[[i]])  <- names(annotations.list)[[i]]
		}

	#4. save in meta.data$annotation
		SO.integrated[["annotation"]] <- Idents(SO.integrated)
		DimPlot(SO.integrated, label=T, label.size=5)

	#5. remove cycling
		SO.integrated <- subset(SO.integrated, idents = "cycling", invert = TRUE)
	
	#6. reorder idents and create DimPlots 
		Idents(SO.integrated) <- factor(Idents(SO.integrated), levels = c("IgA plasma cell","IgG plasma cell","IgM plasma cell","IgD plasma cell","B cell","hepatocyte","cholangiocyte","endothelial cell","stromal cell","T cell","NK CD56bright","NK CD56dim","myeloid cell"))
	
	#7. save 
		# saveRDS(SO.integrated, file="saveRDS/2liver_Ramachandran_Nkongolo.Rds")      #save
	
	#8. plot
		#colors
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
					color.vector <- c(color.vector[1:6], color.vector[11], color.vector[c(7:10,12:13)])
			}
			color.vector.alpha2 <- alpha(color.vector,0.2)	
			color.vector.alpha5 <- alpha(color.vector,0.5)	 
 
	#UMAP representation plots
		p=DimPlot(SO.integrated, cols= color.vector, pt.size=1)
		p
		if(save.plots==TRUE){ggsave(p, file= "output/plots/2liver_DimPlot_legend.tiff", width = 3, height=4)}

		p=DimPlot(SO.integrated, cols= color.vector, pt.size=0.5) + theme_void() + theme(legend.position = 'none')
		p
		if(save.plots==TRUE){ggsave(p, file= "output/plots/2liver_DimPlot_void_pt0-5.tiff", width = 4, height=4)}
		
		p=DimPlot(SO.integrated, cols= color.vector.alpha2, pt.size=0.5) + theme_void() + theme(legend.position = 'none')
		p
		if(save.plots==TRUE){  ggsave(p, file= "output/plots/2liver_DimPlot_void_a2.tiff", width = 4, height=4)  }

	#DotPlots
		dotplot.markers <- c("MZB1","TNFRSF17","MS4A1","ALB","TTR","SOX9","KRT19","CD34","PECAM1","MYL9","COL1A1","TAGLN","COL3A1","CD2","CD3D","NCAM1","XCL1","XCL2","CD160","GNLY","GZMB","FCGR3A","FCN1","VCAN","CLEC10A","MARCO")
		p=DotPlot(SO.integrated, cols = "RdYlBu",features = dotplot.markers) + RotatedAxis() & theme(axis.title.x = element_blank(),axis.title.y = element_blank())
		p
		if(save.plots==TRUE){ggsave(p, file= "output/plots/2liver_dotplot_sc_markers.tiff", width = 10, height=4) }
		
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
		if(save.plots==TRUE){ggsave(p, file= "output/plots/2liver_clustertable.tiff", width = 15, height=7) }
		
		
	# Zoom in on coordinates of plasma cells
		SO.zoom <- subset(SO.integrated, idents=c("IgA plasma cell","IgM plasma cell","IgG plasma cell","IgD plasma cell","IgD plasma cell"))

		#cut everything beyond these UMAP coordinates
			UMAP.coordinates.top    <- -1
			UMAP.coordinates.bottom <- -5
			UMAP.coordinates.left   <- 3
			UMAP.coordinates.right  <- 6
			
			remove.cell.top    <- which(SO.zoom[["umap"]]@cell.embeddings[,2] > UMAP.coordinates.top) %>% names
			remove.cell.bottom <- which(SO.zoom[["umap"]]@cell.embeddings[,2] < UMAP.coordinates.bottom) %>% names
			remove.cell.left   <- which(SO.zoom[["umap"]]@cell.embeddings[,1] < UMAP.coordinates.left) %>% names
			remove.cell.right  <- which(SO.zoom[["umap"]]@cell.embeddings[,1] > UMAP.coordinates.right) %>% names
			
			SO.zoom <- subset(SO.zoom, cells= unique(c(remove.cell.top,remove.cell.bottom,remove.cell.left,remove.cell.right)), invert=TRUE)

			DimPlot(SO.zoom, cols=color.vector[1:4])
		
		p=DimPlot(SO.zoom, cols=color.vector.alpha5[1:4], pt.size=1) +theme_void() + NoLegend()
		p
		if(save.plots==TRUE){  ggsave(p, file= "output/plots/2liver_DimPlot_Ig_zoom.tiff", width = 2, height=2)  }

		#visualise separately
			color.vector.alpha8 <- alpha(color.vector,0.8)	 

			clust.ident=("IgA plasma cell")
			p=DimPlot(SO.zoom, cols.highlight=color.vector.alpha8[1], pt.size=1, sizes.highlight=1, cells.highlight = colnames(subset(SO.zoom,ident=clust.ident))) +theme_void() + NoLegend() 
			if(save.plots==TRUE){  ggsave(p, file= "output/plots/2liver_DimPlot_IgA_zoom.tiff", width = 2, height=2)  }
			
			clust.ident=("IgG plasma cell")
			p=DimPlot(SO.zoom, cols.highlight=color.vector.alpha8[2], pt.size=1, sizes.highlight=1, cells.highlight = colnames(subset(SO.zoom,ident=clust.ident))) +theme_void() + NoLegend()
			if(save.plots==TRUE){  ggsave(p, file= "output/plots/2liver_DimPlot_IgG_zoom.tiff", width = 2, height=2)  }
			
			clust.ident=("IgM plasma cell")
			p=DimPlot(SO.zoom, cols.highlight=color.vector.alpha8[3], pt.size=1, sizes.highlight=1, cells.highlight = colnames(subset(SO.zoom,ident=clust.ident))) +theme_void() + NoLegend()
			if(save.plots==TRUE){  ggsave(p, file= "output/plots/2liver_DimPlot_IgM_zoom.tiff", width = 2, height=2)  }
			
			clust.ident=("IgD plasma cell")
			p=DimPlot(SO.zoom, cols.highlight=color.vector.alpha8[4], pt.size=1, sizes.highlight=1, cells.highlight = colnames(subset(SO.zoom,ident=clust.ident))) +theme_void() + NoLegend()
			if(save.plots==TRUE){  ggsave(p, file= "output/plots/2liver_DimPlot_IgD_zoom.tiff", width = 2, height=2)  }


	# create figure dotplot B and plasma cell markers
	library(RColorBrewer)
		SO.subset <- subset(SO.integrated, idents=c("IgA plasma cell","IgM plasma cell","IgG plasma cell","IgD plasma cell"))
		Idents(SO.subset) <- factor(Idents(SO.subset), levels = c("IgD plasma cell","IgM plasma cell","IgG plasma cell","IgA plasma cell"))
	
		dotplot.features <- c("IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","IGHM","IGHD")
		
		p=DotPlot(SO.subset, cols = "RdYlBu",features = dotplot.features) + RotatedAxis() &  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
		p




# sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=Dutch_Belgium.1252  LC_CTYPE=Dutch_Belgium.1252    LC_MONETARY=Dutch_Belgium.1252 LC_NUMERIC=C                   LC_TIME=Dutch_Belgium.1252    
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] dplyr_1.1.1        patchwork_1.1.1    RColorBrewer_1.1-2 ggplot2_3.4.4      ggsci_2.9          harmony_0.1.0      Rcpp_1.0.7         SeuratObject_4.0.4 Seurat_4.0.2      
# 
# loaded via a namespace (and not attached):
#  [1] nlme_3.1-149          matrixStats_0.58.0    spatstat.sparse_3.0-2 RcppAnnoy_0.0.18      httr_1.4.2            sctransform_0.3.2     tools_4.0.3          
#  [8] utf8_1.2.1            R6_2.5.0              irlba_2.3.3           rpart_4.1-15          KernSmooth_2.23-17    uwot_0.1.10           mgcv_1.8-33          
# [15] lazyeval_0.2.2        colorspace_2.0-0      withr_2.5.0           tidyselect_1.2.0      gridExtra_2.3         compiler_4.0.3        cli_3.6.2            
# [22] plotly_4.9.3          labeling_0.4.2        scales_1.2.1          lmtest_0.9-38         spatstat.data_3.0-1   ggridges_0.5.3        pbapply_1.4-3        
# [29] goftest_1.2-2         stringr_1.5.1         digest_0.6.27         spatstat.utils_3.0-3  pkgconfig_2.0.3       htmltools_0.5.2       parallelly_1.32.0    
# [36] fastmap_1.1.0         htmlwidgets_1.5.3     rlang_1.1.3           shiny_1.6.0           farver_2.1.0          generics_0.1.3        zoo_1.8-12           
# [43] jsonlite_1.7.2        ica_1.0-2             magrittr_2.0.1        Matrix_1.3-4          munsell_0.5.0         fansi_0.4.2           abind_1.4-5          
# [50] reticulate_1.18       lifecycle_1.0.3       stringi_1.5.3         MASS_7.3-58.3         Rtsne_0.15            plyr_1.8.6            grid_4.0.3           
# [57] parallel_4.0.3        listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1         miniUI_0.1.1.1        deldir_1.0-6          lattice_0.20-41      
# [64] cowplot_1.1.1         splines_4.0.3         tensor_1.5            pillar_1.9.0          igraph_1.2.6          spatstat.geom_3.2-1   future.apply_1.7.0   
# [71] reshape2_1.4.4        codetools_0.2-16      leiden_0.3.7          glue_1.6.2            data.table_1.14.0     png_0.1-7             vctrs_0.6.5          
# [78] httpuv_1.5.5          gtable_0.3.0          RANN_2.6.1            purrr_1.0.2           spatstat.core_2.0-0   polyclip_1.10-0       tidyr_1.3.1          
# [85] scattermore_0.7       future_1.26.1         mime_0.10             xtable_1.8-4          RSpectra_0.16-0       later_1.1.0.1         survival_3.2-7       
# [92] viridisLite_0.4.2     tibble_3.2.1          cluster_2.1.0         globals_0.15.1        fitdistrplus_1.1-3    ellipsis_0.3.2        ROCR_1.0-11          
# 
	