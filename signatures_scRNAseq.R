###############################################################
#
# plot single cell module scores
#
###############################################################

#note: the .Rds file with the seurat object of 5 integrated human liver samples can be downloaded from https://zenodo.org/records/10149895 ; DOI 10.5281/zenodo.10149894
#note: if you use this data, please cite the original publication: MacParland et al. 2018 Nature Communications, https://doi.org/10.1038/s41467-018-06318-7

#set working directory
	setwd("C:/my_dir")
	
#load libraries
	library(Seurat)
	library(ggplot2)
	library(RColorBrewer)
	library(ggsci)
	library(patchwork)

#load data 
	SO.integrated <- readRDS(file="data_files/scRNAseq_data/liver_scRNAseq.Rds")

#calculate single cell scores
	#1. Set metric
	scsign.metric <- "modulescore" # "modulescore", "RNA",
	search.setting = FALSE #for modulescore only; set TRUE if genes need to be changed into HUGO symbols
	
	#2. Load gmt file with genesets
		genesets.temp <- read.table(file="data_files/gmt_files/signatures.gmt", sep="\t", fill=T)
		genesets <- list()
			for(i in 1:nrow(genesets.temp)){
				genesets[[i]] <- unname(unlist(genesets.temp[i,which(genesets.temp[i,] !="")][-c(1,2)]) )
				names(genesets)[i] <- genesets.temp[i,1]
			}
		rm(genesets.temp)
	
	#3. algorithm
	{
		if(scsign.metric == "RNA"){ #logNormalized UMI counts 
			DefaultAssay(SO.integrated) <- "RNA"
			SCSIGN.matrix <- matrix(0,length(genesets),ncol(SO.integrated))
			rownames(SCSIGN.matrix) <- names(genesets)
			colnames(SCSIGN.matrix) <- colnames(SO.integrated)
			
			feat.vector <- rownames(SO.integrated)
			denominator <- as.numeric(apply(SO.integrated[["RNA"]][1:nrow(SO.integrated),1:ncol(SO.integrated)], 2, sum)) #sum of all UMIs per cell
			
			for(i in 1:nrow(SCSIGN.matrix) ){
				if( is.element(i,seq(1,100000,100))){ print(paste0("calculating geneset ",i," of ",length(genesets)))}
				row.selection <- as.numeric(na.omit(as.numeric(sapply(genesets[[i]], function(x) which(feat.vector == x))))) #which rows are the genes of geneset i
				temp.matrix <- SO.integrated[["RNA"]][row.selection,] # create matrix with features of geneset i for each cell
				SCSIGN.matrix[i,]  <- as.numeric(100*apply(temp.matrix, 2, sum)) / denominator #(sum of UMIs of genes of geneset i in cell j) divided by (sum of all UMIs in cell j)
			}
			SO.integrated[["SCSIGN"]] <- CreateAssayObject(counts =  SCSIGN.matrix)
			DefaultAssay(SO.integrated) <- "SCSIGN"
		}
		
		if(scsign.metric == "modulescore"){ #use AddModuleScore seurat functionality to calculate gene set score
			DefaultAssay(SO.integrated) <- "RNA"
			
			module.feats <- genesets
			
			ncol.meta.data <- length(SO.integrated@meta.data)
			SO.integrated <- AddModuleScore(object = SO.integrated,  features = module.feats,  name = 'module',  assay= 'RNA',  search = search.setting)
			ncol.meta.data.module <- length(SO.integrated@meta.data)
			
				if((ncol.meta.data.module-ncol.meta.data)==length(genesets)){ print("all geneset scores succesfully calculated")}
				if((ncol.meta.data.module-ncol.meta.data)!=length(genesets)){ print("failure: not all geneset scores succesfully calculated")}
		
			SCSIGN.matrix <- t(SO.integrated@meta.data[,-(1:ncol.meta.data)])	
			rownames(SCSIGN.matrix) <- names(genesets)
			
			delete.cols <- names(SO.integrated@meta.data[,-(1:ncol.meta.data)])
			for(i in delete.cols){ #delete module scores from meta.data
				SO.integrated[[i]] <- NULL
			}
			
			#scale SCSIGN.matrix from 0 to 1
				scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
				SCSIGN.matrix <-  t(apply(SCSIGN.matrix,1,scale01))
						
			SO.integrated[["SCSIGN"]] <- CreateAssayObject(counts =  SCSIGN.matrix)
			DefaultAssay(SO.integrated) <- "SCSIGN"
		}
	}
	
	#4. Visualise
		#set colors
			color.vector <- pal_npg("nrc")(10) #requires ggsci package
			color.vector <- c(color.vector,"brown","khaki","orange3","orchid2")
	
		#high-risk genes
			feat = rownames(SO.integrated)[1]
			
			#calculate quantiles
				feat.values <- as.numeric(SO.integrated[[DefaultAssay(SO.integrated)]][rownames(SO.integrated) %in% feat,])
			
				quant10 <- quantile(feat.values,0.10)
				quant25 <- quantile(feat.values,0.25)
				quant75 <- quantile(feat.values,0.75)
				quant90 <- quantile(feat.values,0.90)
			
			#plot 
				p=VlnPlot(SO.integrated, feat, cols= color.vector, pt.size=0.2) + 
					#stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)+
					guides(fill = guide_legend(override.aes = list(size=4), ncol=1)) +
					geom_hline(yintercept = c(quant10, quant90), linetype='dashed', col = 'red', linewidth=1) +
					theme(axis.title.x=element_blank()) +
					theme(legend.position = 'none') +
					ggtitle("")
				p

		#low-risk genes
			feat = rownames(SO.integrated)[2]
			
			#calculate quantiles
				feat.values <- as.numeric(SO.integrated[[DefaultAssay(SO.integrated)]][rownames(SO.integrated) %in% feat,])
			
				quant10 <- quantile(feat.values,0.10)
				quant25 <- quantile(feat.values,0.25)
				quant75 <- quantile(feat.values,0.75)
				quant90 <- quantile(feat.values,0.90)
			
			#plot
				p=VlnPlot(SO.integrated, feat, cols= color.vector, pt.size=0.2) + 
					#stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)+
					guides(fill = guide_legend(override.aes = list(size=4), ncol=1)) +
					geom_hline(yintercept = c(quant10, quant90), linetype='dashed', col = 'red', linewidth=1) +
					theme(axis.title.x=element_blank()) +
					theme(legend.position = 'none') +
					ggtitle("")
				p
			