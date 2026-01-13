#######################################################
#
# Creates heatmaps to accompany the NTP stratification
# showing expression of the 557 gene signature
# This requires the output of the NTP algorithm
#  
########################################################

#set working directory
	setwd("C:/my_dir")
	
#load libraries
	library(circlize)
	library(ComplexHeatmap)

#0 Load gene signature
############################################
	gene.signature <- read.delim(file="output/training_RPM_filtered_signature_result_ordered.txt",header=T,check.names=F)
	risk.genes <- rep(0, nrow(gene.signature))
	risk.genes[which(gene.signature[,2] > 1)] <- "high-risk genes" 
	risk.genes[which(gene.signature[,2] < 1)] <- "low-risk genes"		
	gene.signature$risk <- risk.genes
		
#1 Training set heatmap
############################################
	#load data
		#extract gene expression from TPM data
			exp.dataset <- read.delim(file="output/training_NTP_sorted_dataset.gct",header=T,skip=2,check.names=F)
			rownames(exp.dataset) <- exp.dataset[,1]
			mat <- exp.dataset[,-c(1:2)]
			
		#read meta data
			meta <- read.delim(file="data_files/meta_data/training_full_meta.txt",header=T,check.names=F)

		#Standardize data
			cal_z_score <- function(x){
				(x - mean(x)) / sd(x)
			}
		
			mat.norm <- t(apply(mat, 1, cal_z_score))
	
	#change colors	
		#Users should always use circlize::colorRamp2() function to generate the color mapping function in Heatmap().
		col_fun = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick"))
	
	#add annotations
		#column annotation
			#read cls & change 1/2/3 into poor/intermediate/good
  			cls.ordered <- read.table(file="output/training_NTP_PoorIntGood_predicted_sorted.cls",header=F,check.names=F, skip=2)
  			cls.ordered[cls.ordered == 1] <- "poor"
  			cls.ordered[cls.ordered == 2] <- "good"
  			cls.ordered[cls.ordered == 3] <- "intermediate"
  			cls.ordered <- unlist(cls.ordered)

			# Map meta data to the columns of mat.norm using SubjectName
  			treatment <- meta$Treated[match(colnames(mat.norm), meta$SubjectName)]
  			treatment <- factor(treatment, levels=c(0, 1), labels=c("no", "yes"))
  			
  			HCC <- meta$HCC_event[match(colnames(mat.norm), meta$SubjectName)]
  			HCC <- factor(HCC, levels=c(0, 1), labels=c("no", "yes"))
  			
  			#create column annotations for complexheatmap
    			# Combine column annotations
      			column_ha <- HeatmapAnnotation(
      			  prognosis=cls.ordered,
      			  treated=treatment,
      			  HCC=HCC,
      			  col=list(
      			    prognosis=c("poor"="red", "intermediate"="lightgrey", "good"="blue"),
      			    treated=c("no"="snow1", "yes"="purple"),
      			    HCC=c("no"="snow1", "yes"="salmon")
      			  ),
      			  simple_anno_size=unit(0.3, "cm"),
      			  show_annotation_name=c(prognosis=FALSE, treatment=FALSE, HCC=FALSE)
      			)
    			
			
		#row annotation
			#filter signature for genes that are NOT in the heatmap's dataset
			gene.signature.filtered <- gene.signature[gene.signature$probeid %in% rownames(mat.norm),]
	
			#create row annotations for complexheatmap
			row_ha_left = rowAnnotation(HCC_risk = gene.signature.filtered$risk, 
				col = list(HCC_risk = c("high-risk genes" = "orange", "low-risk genes" = "seagreen4")),
				simple_anno_size = unit(0.3, "cm"),
				show_annotation_name = c(HCC_risk = FALSE),
				annotation_legend_param = list(HCC_risk = list(title = "HCC incidence"))
			)
			
			highlight.genes <- c("IGHA1","IGHA2")
			highlight.genes.row <- match(highlight.genes, rownames(mat.norm))
			row_ha_right = rowAnnotation(highlighted_genes = anno_mark(at = highlight.genes.row, labels = highlight.genes))

	#visualize
		ht=Heatmap(mat.norm, 
			cluster_rows=FALSE, 
			cluster_columns=FALSE, 
			col = col_fun, 
			show_row_names=FALSE, 
			show_column_names=FALSE, 
			top_annotation=column_ha,
			left_annotation=row_ha_left,
			right_annotation=row_ha_right,
			heatmap_legend_param = list(title = "Z-score")
		)
		ht			

		#save pdf
		  pdf("plots_pdf/Fig_B_heatmap_training.pdf", width = 5.5, height = 4)
		  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
		  dev.off()
		  
		  
#2 Validation set heatmap
############################################
	#load data
		#extract gene expression from TPM data
			exp.dataset <- read.delim(file="output/validation_NTP_sorted_dataset.gct",header=T,skip=2,check.names=F)
			rownames(exp.dataset) <- exp.dataset[,1]
			mat <- exp.dataset[,-c(1:2)]

		#read meta data
		  meta <- read.delim(file="data_files/meta_data/validation_full_meta.txt",header=T,check.names=F)
			
		#Standardize data
			cal_z_score <- function(x){
				(x - mean(x)) / sd(x)
			}
		
			mat.norm <- t(apply(mat, 1, cal_z_score))
	
	#change colors	
		#Users should always use circlize::colorRamp2() function to generate the color mapping function in Heatmap().
		col_fun = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick"))
	
	#add annotations
		#column annotation
			#read cls & change 1/2/3 into poor/intermediate/good
			cls.ordered <- read.table(file="output/validation_NTP_PoorIntGood_predicted_sorted.cls",header=F,check.names=F, skip=2)
			cls.ordered[cls.ordered == 1] <- "poor"
			cls.ordered[cls.ordered == 2] <- "good"
			cls.ordered[cls.ordered == 3] <- "intermediate"
			cls.ordered <- unlist(cls.ordered)
			
		# Map meta data to the columns of mat.norm using SubjectName
			treatment <- meta$Treated[match(colnames(mat.norm), meta$SubjectName)]
			treatment <- factor(treatment, levels=c(0, 1), labels=c("no", "yes"))
			
			HCC <- meta$HCC_event[match(colnames(mat.norm), meta$SubjectName)]
			HCC <- factor(HCC, levels=c(0, 1), labels=c("no", "yes"))
			
		#create column annotations for complexheatmap
			# Combine column annotations
			column_ha <- HeatmapAnnotation(
			  prognosis=cls.ordered,
			  treated=treatment,
			  HCC=HCC,
			  col=list(
			    prognosis=c("poor"="red", "intermediate"="lightgrey", "good"="blue"),
			    treated=c("no"="snow1", "yes"="purple"),
			    HCC=c("no"="snow1", "yes"="salmon")
			  ),
			  simple_anno_size=unit(0.3, "cm"),
			  show_annotation_name=c(prognosis=FALSE, treatment=FALSE, HCC=FALSE)
			)
			
			
		#row annotation
			#filter signature for genes that are NOT in the heatmap's dataset
			gene.signature.filtered <- gene.signature[gene.signature$probeid %in% rownames(mat.norm),]
	
				#create row annotations for complexheatmap
				row_ha_left = rowAnnotation(HCC_risk = gene.signature.filtered$risk, 
					col = list(HCC_risk = c("high-risk genes" = "orange", "low-risk genes" = "seagreen4")),
					simple_anno_size = unit(0.3, "cm"),
					show_annotation_name = c(HCC_risk = FALSE),
					annotation_legend_param = list(HCC_risk = list(title = "HCC incidence")) 
				)
			
			highlight.genes <- c("IGHA1","IGHA2")
			highlight.genes.row <- match(highlight.genes, rownames(mat.norm))
			row_ha_right = rowAnnotation(highlighted_genes = anno_mark(at = highlight.genes.row, labels = highlight.genes))
		
	#visualize
		ht=Heatmap(mat.norm, 
			cluster_rows=FALSE, 
			cluster_columns=FALSE, 
			col = col_fun, 
			show_row_names=FALSE, 
			show_column_names=FALSE, 
			top_annotation=column_ha,
			left_annotation=row_ha_left,
			right_annotation=row_ha_right,
			heatmap_legend_param = list(title = "Z-score")
		)
		ht			
			
		#save pdf
		pdf("plots_pdf/Fig_2E_heatmap_validation.pdf", width = 5.5, height = 4)
		draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
		dev.off()


