############################################
#
#	create IgA IHC staining heatmap
#
############################################

#set working directory
	setwd("C:/my_dir")

#load libraries
	library(ComplexHeatmap)
	library(Seurat)
	library(ggplot2)
	library(RColorBrewer)
	library(ggsci)
	library(patchwork)
	
#load data sets
	IHC.df <- read.table(file="data_files/IHC_staining_data/IgA_staining.txt", header=T, sep="\t")
	exp.dataset <- read.delim(file="data_files/collapsed_data/GSE237330_training_RPM.gct",header=T,skip=2,check.names=F)
		rownames(exp.dataset) <- exp.dataset[,1]
		exp.df <- exp.dataset[,-c(1:2)]

#log10 transform
	exp.log10.df <- log(1+exp.df,10)

#calculate marker genes
	#load data
		#note: the .Rds file with the seurat object of 5 integrated human liver samples can be downloaded from https://zenodo.org/records/10149895 ; DOI 10.5281/zenodo.10149894
		#note: is you use this data, please cite the original publication: MacParland et al. 2018 Nature Communications, https://doi.org/10.1038/s41467-018-06318-7

		SO.integrated <- readRDS(file="data_files/scRNAseq_data/liver_scRNAseq.Rds")
	
	#calculate markers for IgA plasma cells
		min.pct.expression = 0.25 #standard setting: 0.25
		min.logfc = 0.25 #0.25 is standard
		p.val.cutoff <- (1/10^3) #eg. (1/10^3) 
		
		DefaultAssay(SO.integrated) <- "RNA"
			IgA.markers.df <- FindMarkers(SO.integrated, ident.1 = "IgA plasma cell", min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
			IgA.markers.df <- IgA.markers.df[which(IgA.markers.df$p_val_adj < (p.val.cutoff) ),] #select genes with p_val_adj > p.val.cutoff setting
		
		IgA.markers <- rownames(IgA.markers.df)
		DotPlot(SO.integrated, cols = "RdYlBu",features = rev(IgA.markers)) + RotatedAxis() 


		IgA.markers <- IgA.markers[1:20]
		#IgA.markers <- c("AL928768.3","IGHA2","MZB1","DERL3","TNFRSF17","JCHAIN","CPNE5","POU2AF1","IGHA1", "FCRL5", "JSRP1", "IGLC3", "ABCB9","RP11-16E12.2","CD79A","ITM2C","PNOC","IGLV6-57","SPAG4","FKBP11","CD27","PIM2")       

	my.markers <- IgA.markers[IgA.markers %in% rownames(exp.df)]
	exp.df <- exp.df[rownames(exp.df) %in% my.markers,] #filter markers for those that are not in the gene expression database
	exp.log10.df <- exp.log10.df[rownames(exp.log10.df) %in% my.markers,] #filter markers for those that are not in the gene expression database


#order gene expression matrix
	#order X-axis according to IgA pos/mm²
		posmm2.df <- IHC.df[, colnames(IHC.df) %in% c("SubjectName","Posmm2")]
		posmm2.df <- posmm2.df[ order(posmm2.df$Posmm2),]
		exp.df <- exp.df[, match(posmm2.df$SubjectName, colnames(exp.df))]

	#order X-axis according to order of marker genes
		exp.df <- exp.df[match(my.markers, rownames(exp.df)),]

#Standardize data
	cal_z_score <- function(x){
		(x - mean(x)) / sd(x)
	}
	exp.df.norm <- t(apply(exp.df, 1, cal_z_score))

#calculate correlations and pvalues for all genes (log10 transformed)
	#select IHC samples
		exp.log10.df <- exp.log10.df[,colnames(exp.log10.df) %in% IHC.df$SubjectName]
		exp.log10.df <- exp.log10.df[,match(IHC.df$SubjectName, colnames(exp.log10.df))]
		
		if(is.element(FALSE,colnames(exp.log10.df) == IHC.df$SubjectName)){print("error: IHC.df$SubjectName does not match expression matrix colnames")} #QC

		cor.matrix <- matrix(0,nrow(exp.log10.df),3)
		colnames(cor.matrix) <- c("R","pvalue","slope")
		rownames(cor.matrix) <- rownames(exp.log10.df)
		
		for(i in 1:nrow(exp.log10.df)){
			cor.matrix[i,1]	<- round(cor(IHC.df$Posmm2, as.numeric(exp.log10.df[i,])),2) #correlation value
			cor.matrix[i,2]	<- cor.test(IHC.df$Posmm2, as.numeric(exp.log10.df[i,]))$p.value #pvalue
			cor.matrix[i,3]	<- as.numeric(lm(IHC.df$Posmm2 ~ as.numeric(exp.log10.df[i,]))$coefficient[2]) #slope, this is change in Y (continuous meta data factor) per unit increase in X (gene expression).
		}
		cor.matrix.pos <- cor.matrix[cor.matrix[,1]>=0,]
		cor.matrix.neg <- cor.matrix[cor.matrix[,1]<0,]
		
		cor.matrix.pos <- cor.matrix.pos[order(cor.matrix.pos[,1], decreasing=TRUE),]
		cor.matrix.neg <- cor.matrix.neg[order(cor.matrix.neg[,1], decreasing=TRUE),]
		
		cor.matrix.ordered <- cor.matrix[match(rownames(exp.df),rownames(cor.matrix)),]
		
		#change pvalues to *
			p.values.IgAmarkers <- as.numeric(cor.matrix.ordered[,colnames(cor.matrix.ordered) %in% "pvalue"])
			p.star <- as.numeric(cor.matrix.ordered[,colnames(cor.matrix.ordered) %in% "pvalue"])

			p.star[p.star < 0.001] <- "***"
			p.star[p.star < 0.05 & p.star > 0.01]  <- "*"
			p.star[p.star < 0.01 & p.star > 0.001]  <- "**"
			p.star[p.star >= 0.05] <- ""
		
			cor.matrix.ordered <- cbind(cor.matrix.ordered,p.star)
			colnames(cor.matrix.ordered)[4] <- "pstar"
			cor.matrix.ordered

#specific genes correlation plots
	if(is.element(FALSE,colnames(exp.log10.df) == IHC.df$SubjectName)){print("error: IHC.df$SubjectName does not match expression matrix colnames")} #QC

	#plot CD79A
		IHC.df$CD79A <- as.numeric(unname(exp.log10.df[rownames(exp.log10.df) %in% "CD79A",]))

		cor.val <- round(cor(IHC.df$Posmm2, IHC.df$CD79A), 2)
		cor.p <- signif(cor.test(IHC.df$Posmm2, IHC.df$CD79A)$p.value,2)
		slope <- lm(IHC.df$Posmm2 ~ IHC.df$CD79A)

		p=ggplot(IHC.df, aes(x = CD79A, y = Posmm2)) + 
			geom_point(size=3) +
			theme_classic() +
			geom_smooth(method='lm', se=FALSE, color="black") +
			annotate(x = 0, y = 8,   hjust=0, geom = "text", label = paste0("Pearson's r=",cor.val), size = 6)+
			annotate(x = 0, y = 7.5, hjust=0, geom = "text", label = paste0("p-value=",cor.p), size = 6) +
			xlab("CD79A (log10 RPM)") + ylab("IgA+ cells / mm²") +
			theme(axis.text.x = element_text(size = 20, color="black")) +
			theme(axis.text.y = element_text(size = 20, color="black")) +
			theme(axis.title.x = element_text(size = 20)) +
			theme(axis.title.y = element_text(size = 20)) +
			theme(legend.title = element_blank()) +
			theme(legend.text = element_text(size = 20)) 
		p

	#plot IGHA2
		IHC.df$IGHA2 <- as.numeric(unname(exp.log10.df[rownames(exp.log10.df) %in% "IGHA2",]))

		cor.val <- round(cor(IHC.df$Posmm2, IHC.df$IGHA2), 2)
		cor.p <- signif(cor.test(IHC.df$Posmm2, IHC.df$IGHA2)$p.value,2)
		slope <- lm(IHC.df$Posmm2 ~ IHC.df$IGHA2)

		p=ggplot(IHC.df, aes(x = IGHA2, y = Posmm2)) + 
			geom_point(size=3) +
			theme_classic() +
			geom_smooth(method='lm', se=FALSE, color="black") +
			annotate(x = 0, y = 8,   hjust=0, geom = "text", label = paste0("Pearson's r=",cor.val), size = 6)+
			annotate(x = 0, y = 7.5, hjust=0, geom = "text", label = paste0("p-value=",cor.p),   size = 6)+
			xlab("IGHA2 (log10 RPM)") + ylab("IgA+ cells / mm²") +
			theme(axis.text.x = element_text(size = 20, color="black")) +
			theme(axis.text.y = element_text(size = 20, color="black")) +
			theme(axis.title.x = element_text(size = 20)) +
			theme(axis.title.y = element_text(size = 20)) +
			theme(legend.title = element_blank()) +
			theme(legend.text = element_text(size = 20)) 
		p

#plot POU2AF1
		IHC.df$POU2AF1 <- as.numeric(unname(exp.log10.df[rownames(exp.log10.df) %in% "POU2AF1",]))

		cor.val <- round(cor(IHC.df$Posmm2, IHC.df$POU2AF1), 2)
		cor.p <- signif(cor.test(IHC.df$Posmm2, IHC.df$POU2AF1)$p.value,2)
		slope <- lm(IHC.df$Posmm2 ~ IHC.df$POU2AF1)

		p=ggplot(IHC.df, aes(x = POU2AF1, y = Posmm2)) + 
			geom_point(size=3) +
			theme_classic() +
			geom_smooth(method='lm', se=FALSE, color="black") +
			annotate(x = 0, y = 8,   hjust=0, geom = "text", label = paste0("Pearson's r=",cor.val), size = 6)+
			annotate(x = 0, y = 7.5, hjust=0, geom = "text", label = paste0("p-value=",cor.p), size = 6)+
			xlab("POU2AF1 (log10 RPM)") + ylab("IgA+ cells / mm²") +
			theme(axis.text.x = element_text(size = 20, color="black")) +
			theme(axis.text.y = element_text(size = 20, color="black")) +
			theme(axis.title.x = element_text(size = 20)) +
			theme(axis.title.y = element_text(size = 20)) +
			theme(legend.title = element_blank()) +
			theme(legend.text = element_text(size = 20)) 
		p
		
#visualize in heatmap
	#change colors	of heatmap
		#Users should always use circlize::colorRamp2() function to generate the color mapping function in Heatmap().
		library(circlize)
		col_fun = colorRamp2(c(-3, 0, 3), c("navy", "white", "firebrick"))
	
	#add annotations		
		#column annotation: IHC IgA+ cells/mm²		
			posmm2.vector <- posmm2.df[,2]
			IHC.color = colorRamp2(c(0, max(round(posmm2.df[,2]))), c("white", "purple"))
		
		#column annotation: Metavir
			metavir.vector 	<- IHC.df$Metavir[ match( colnames(exp.df.norm), IHC.df$SubjectName) ]
			metavir.color = colorRamp2(c(0,1,2,3,4), c("white", "snow1","khaki","orange3","brown"))		
		
		#add all column annotations
			column_ha = HeatmapAnnotation(
					"Metavir" = metavir.vector,
					"IgA+ cells / mm²" = posmm2.vector,
					col = list("IgA+ cells / mm²" = IHC.color, "Metavir" = metavir.color),
					simple_anno_size = unit(0.5, "cm")
			)
		

		#add pvalue right annotations
			#for numbers
				#right_ha = rowAnnotation(
					#"p-value" = anno_text(signif(cor.matrix.ordered[,colnames(cor.matrix.ordered) %in% "pvalue"],2))
				#)
			
			#for * / ** / ***			
				right_ha = rowAnnotation(
					"p-value" = anno_text(cor.matrix.ordered[,colnames(cor.matrix.ordered) %in% "pstar"])
				)
			
			#reproduce heatmap with metavir and IgA+ cell stain count
				ht=Heatmap(exp.df.norm,
					cluster_rows=FALSE,
					cluster_columns=FALSE,
					col = col_fun,
					show_row_names=TRUE,
					show_column_names=TRUE,
					right_annotation = right_ha,
					top_annotation=column_ha,
					heatmap_legend_param = list(title = "Z-score")
					)
				ht	
		
		