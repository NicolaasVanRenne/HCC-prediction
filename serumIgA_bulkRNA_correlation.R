###########################################################
#
# correlate serum IgA with bulk RNAseq
# and project on liver scRNA-seq 
#
##########################################################

#set working directory
	setwd("C:/my_dir")

#load libraries
	library(ggplot2)
	library(dplyr)
	library(ggrepel)
	library(Seurat)

#correlate bulkseq RNA with serum IgA
#############################################

#load databases
	IgA.df 	<- read.table(file="data_files/serum_data/serumIgA_bulkRNA.txt", sep="\t", header=T) 
	
	#format SerumIgA as numeric
		IgA.df$SerumIgA <- as.numeric(IgA.df$SerumIgA)
		
	#log transform
		IgA.df$IGHA1 <- log10(1+IgA.df$IGHA1)
		IgA.df$IGHA2 <- log10(1+IgA.df$IGHA2)
		
#calculate correlations between serum IgA and all bulkseq genes
	#1. import training_56 data
		my.file = "data_files/collapsed_data/GSE237330_training_RPM.gct"
		read.gct<-function(filename="NULL"){
			if (regexpr(".gct$",filename)==-1){
				stop("### input data should be .gct file! ###")
			}
			data<-read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T,check.names=F)
			data<-data[-1]
			return(data)
		}
			
		exp.df      <- read.gct(file=my.file)
	
	#2. transform log10
		exp.df <- log10(exp.df + 1)
		
	#3. filter for samples that have serum IgA	
		exp.df <- exp.df[,colnames(exp.df) %in% IgA.df$SubjectName]

	#4. reorder Serum IgA dataset and bulk RNAseq dataset
		exp.df <- exp.df[,order(colnames(exp.df))]
		IgA.df <- IgA.df[order(IgA.df$SubjectName),]
	
	#5. calculate correlations and pvalues for all genes
		cor.matrix <- matrix(0,nrow(exp.df),3)
		colnames(cor.matrix) <- c("R","pvalue","slope")
		rownames(cor.matrix) <- rownames(exp.df)
		
		for(i in 1:nrow(exp.df)){
			cor.matrix[i,1]	<- round(cor(IgA.df$SerumIgA, as.numeric(exp.df[i,])),2) #correlation value
			cor.matrix[i,2]	<- cor.test(IgA.df$SerumIgA, as.numeric(exp.df[i,]))$p.value #pvalue
			cor.matrix[i,3]	<- as.numeric(lm(IgA.df$SerumIgA ~ as.numeric(exp.df[i,]))$coefficient[2]) #slope, this is change in Y (continuous meta data factor) per unit increase in X (gene expression).
		}
		cor.matrix.pos <- cor.matrix[cor.matrix[,1]>=0,]
		cor.matrix.neg <- cor.matrix[cor.matrix[,1]<0,]
		
		cor.matrix.pos <- cor.matrix.pos[order(cor.matrix.pos[,1], decreasing=TRUE),]
		cor.matrix.neg <- cor.matrix.neg[order(cor.matrix.neg[,1], decreasing=TRUE),]
		
		cor.matrix.ordered <- rbind(cor.matrix.pos,cor.matrix.neg)
			
#Publication quality volcano plot
	volcano.df <- as.data.frame(cor.matrix.ordered)

	#add column of colored genes (top of two branches of volcano plot)
		volcano.df$Direction <- "None"
		volcano.df$Direction[volcano.df$slope > 0 & volcano.df$pvalue < 0.025] <- "Serum IgA correlated" 
		volcano.df$Direction[volcano.df$slope < 0 & volcano.df$pvalue < 0.025] <- "Serum IgA anti-correlated"
		volcano.df$Direction <- as.factor(volcano.df$Direction)
		volcano.df$Direction <- factor(volcano.df$Direction, levels = c("Serum IgA anti-correlated", "None", "Serum IgA correlated"))

	# partially transparent points by setting `alpha`
		linecolors <- c("navy", "darkgrey", "firebrick")
		fillcolors <- c("navy", "darkgrey", "firebrick")

	#add labels
		volcano.df$genelabel <- NA
		
		DN.genes <- NULL 
		UP.genes <- c("IGHA2")
		
		volcano.df$genelabel[rownames(volcano.df) %in% DN.genes] <- DN.genes
		volcano.df$genelabel[rownames(volcano.df) %in% UP.genes] <- UP.genes


		options(ggrepel.max.overlaps = Inf)
		p=ggplot(volcano.df,  aes(x=slope, y=-log10(pvalue), col=Direction, fill=Direction, label=genelabel)) + 
			geom_point(shape = 21, alpha = 0.2, size = 1) +
			scale_color_manual(values=linecolors) +
			scale_fill_manual(values=fillcolors) +
			geom_hline(yintercept=(log10(0.025)*-1), linetype='dashed', col = 'black') +
			labs(y= "-log10(p-value)", x = "change in Serum IgA per unit \n increase in gene expression") +
			theme_classic() +
			theme(axis.text.x = element_text(size = 20, color="black")) +
			theme(axis.text.y = element_text(size = 20, color="black")) +
			theme(axis.title.x = element_text(size = 20)) +
			theme(axis.title.y = element_text(size = 20)) +
			#theme(legend.title = element_blank()) +
			theme(legend.text = element_text(size = 20)) +
			theme(legend.position="top") +guides(col=guide_legend(override.aes = list(size=5), nrow=3)) +
			geom_text_repel(
				#color		 = "black",
				size = 6,
				segment.size = 0.5, #line thickness
				force        = 1,
				nudge_x      = 500,
				direction    = "y",
				hjust        = 1, #horizontal justification:  0 means left-justified, 0.5 means centered, and 1 means right-justified.
				show.legend=FALSE
			) 
		p


# project correlating genes on scRNAseq as module scores 
############################################################

#load scRNAseq dataset of the liver 
	#note: the .Rds file with the seurat object of 5 integrated human liver samples can be downloaded from https://zenodo.org/records/10149895 ; DOI 10.5281/zenodo.10149894
	#note: is you use this data, please cite the original publication: MacParland et al. 2018 Nature Communications, https://doi.org/10.1038/s41467-018-06318-7

	SO.integrated <- readRDS(file="data_files/scRNAseq_data/liver_scRNAseq.Rds")


#calculate module scores
	###1. Settings
		search.setting = FALSE #for modulescore only; set TRUE if genes need to be searched for alt symbols
	
	###2. Set gene set file
		genesets <- list()
		genesets[[1]] <- rownames(cor.matrix.pos[cor.matrix.pos[,2] <0.025,]) 	; names(genesets)[1] <- "PosCorP025"
		genesets[[2]] <- rownames(cor.matrix.neg[cor.matrix.neg[,2] <0.025,]) 	; names(genesets)[2] <- "NegCorP025"					
	
	###3. algorithm and store in [["SCSIGN"]]
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
	
	###4. Visualise
		#create color vector
			library(ggsci)
		
			color.vector <- pal_npg("nrc")(10)
			color.vector <- c(color.vector,"brown","khaki","orange3","orchid2")
			color.vector.alpha2 <- alpha(color.vector,0.2)						
		
		#Serum IgA correlated
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
					ggtitle("Serum IgA correlated")
				p
	
		#Serum IgA anti-correlated
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
					ggtitle("Serum IgA anti-correlated")
				p

