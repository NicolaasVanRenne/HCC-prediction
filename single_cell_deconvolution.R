###############################################################
#
#	Deconvolution of bulk RNAseq using MuSiC algorithm
# 	https://xuranw.github.io/MuSiC/articles/MuSiC.html
#	cite: Wang X. et al. 2019 Nature Communications
#
###############################################################

#set working directory
	setwd("C:/my_dir")
	
# Load packages 
	library(dplyr)
	library(Seurat)
	library(ggplot2)
	library(RColorBrewer)
	library(MuSiC)
	library(Biobase)
	library(ggpubr)
	library(ggsci)

#set memory
	options(future.globals.maxSize = 32000 * 1024^2)
	
#increase memory limit
		memory.limit(999999999999)


# Part 1: load data and gene expression sets
#############################################
	#1. load scRNAseq dataset of the liver 
		#note: the .Rds file with the seurat object of 5 integrated human liver samples can be downloaded from https://zenodo.org/records/10149895 ; DOI 10.5281/zenodo.10149894
		#note: is you use this data, please cite the original publication: MacParland et al. 2018 Nature Communications, https://doi.org/10.1038/s41467-018-06318-7

		SO.integrated <- readRDS(file="data_files/scRNAseq_data/liver_scRNAseq.Rds")
			
	#2. calculate real proportions ('ground truth') for benchmarking pseudobulk estimates 
	{
		real <- table(Idents(SO.integrated), SO.integrated@meta.data$orig.ident) 
		denominator <-  as.numeric(colSums(real))	
		real.normalized <- real
		for(i in 1:ncol(real)){
			real.normalized[,i] <- real[,i]/denominator[i]
		}
		real.normalized <- t(real.normalized)
	}
	
	#3. Create pseudobulk seurat object
	{
		SO.pseudobulk <- SO.integrated
		orig.ident.vector <- unique(SO.integrated$orig.ident)
		
		SO.list <- list()
		for(i in seq_along(orig.ident.vector)){
			SO.list[[i]] <- subset(SO.pseudobulk, subset = orig.ident == orig.ident.vector[i])
		}
		
		bulk.matrix <- c()
		for(i in seq_along(orig.ident.vector)){
			bulk.matrix   <-	cbind(bulk.matrix,rowSums(as.matrix(SO.list[[i]][["RNA"]]@counts[,]))) #RNA raw counts 
		}
		
		bulk.matrix <- bulk.matrix[-which(rowSums(bulk.matrix) == 0),]
		
		colnames(bulk.matrix) <- orig.ident.vector
		
		SO.pseudobulk <- CreateSeuratObject(counts=bulk.matrix, min.cells=3, min.features = 200)
		SO.pseudobulk[["orig.ident"]] <- "Pseudobulk"
	}
	
	#4. Pseudobulk SeuratToExpressionSet 
	{
		individual.ids <-  colnames(SO.pseudobulk) 
		print(paste("there are", length(unique(individual.ids)),"unique individuals:", paste(unique(individual.ids),collapse=" , "))) #print as quality control
		
		base::names(individual.ids) <- colnames(SO.pseudobulk)
		individual.ids <- base::factor(individual.ids)
		n.individuals <- base::length(base::levels(individual.ids))
		sample.ids <- colnames(SO.pseudobulk)
		pseudobulk.pheno <- base::data.frame(check.names=F, check.rows=F,
									stringsAsFactors=F,
									row.names=sample.ids,
									SubjectName=individual.ids,
									cellType=Idents(SO.pseudobulk))
			
		pseudobulk.meta <- base::data.frame(labelDescription=base::c("SubjectName","cellType"),row.names=base::c("SubjectName","cellType"))
		pseudobulk.pdata <- methods::new("AnnotatedDataFrame",data=pseudobulk.pheno,varMetadata=pseudobulk.meta)
		pseudobulk.data <- base::as.matrix(SO.pseudobulk[["RNA"]]@counts[,sample.ids,drop=F])
		pseudobulk.eset <- Biobase::ExpressionSet(assayData=pseudobulk.data,phenoData=pseudobulk.pdata)
	}
	
	#5. Single Cell SeuratToExpressionSet 
	{
		individual.ids <- SO.integrated$orig.ident 
		print(paste("there are", length(unique(individual.ids)),"unique individuals:", paste(unique(individual.ids),collapse=" , "))) #print as quality control
		
		base::names(individual.ids) <- colnames(SO.integrated)
		individual.ids <- base::factor(individual.ids)
		n.individuals <- base::length(base::levels(individual.ids))
		sample.ids <- colnames(SO.integrated)
		sc.pheno <- base::data.frame(check.names=F, check.rows=F,
									stringsAsFactors=F,
									row.names=sample.ids,
									SubjectName=individual.ids,
									cellType=Idents(SO.integrated))
								
		sc.meta <- base::data.frame(labelDescription=base::c("SubjectName","cellType"),row.names=base::c("SubjectName","cellType"))
		sc.pdata <- methods::new("AnnotatedDataFrame",data=sc.pheno,varMetadata=sc.meta)
		sc.eset <- Biobase::ExpressionSet(assayData=base::as.matrix(SO.integrated[["RNA"]]@counts[,sample.ids,drop=F]),phenoData=sc.pdata) 
	}
	

# Part 2: Benchmark MuSiC (pre-grouped) on pseudobulk 
#####################################################

	#1. Manually pre-group clusters based on compartment (immune cells, epithelial, endothelial, stromal)
		clusters.type = list(
			C1  = c("NK CD56dim","NK CD56bright","T cell","kupffer cell","cDC","monocyte","IgA plasma cell","IgG plasma cell","IgM plasma cell","B cell"),
			C2  = c("cholangiocyte","hepatocyte"),
			C3  = c("endothelial cell"), 
			C4  = c("stromal cell")
			)

	#2. run MuSiC algorithm 	
	{	
		# Create clusters object
			cl.type = as.character(sc.eset$cellType)
			for(cl in 1:length(clusters.type)){
				cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
			}
			pData(sc.eset)$clusterType = factor(cl.type, levels = c(names(clusters.type)))
			
			# show number of selected cell types
			unlist(clusters.type)
	
		# Create IEmarkers: a list of markers of within group differentially expressed genes
			{
			# settings
				min.pct.expression = 0.50 #minimal expression of marker genes in single cell celltype
				min.logfc = 0.25 # minimal logfc diffexp of marker genees
				p.val.cutoff <- (1/10^3) # p-value cutoff of diffexp of marker genes
				DefaultAssay(SO.integrated) <- "RNA" 
			
			# run find markers algorithm
				SO.backup <- SO.integrated
				IEmarkers <- vector(mode='list', length=length(clusters.type))
			
				for(j in 1:length(clusters.type)){
					print(paste("calculating within group marker genes for group",j))
					SO.integrated <- subset(SO.backup, idents = clusters.type[[j]])
					cluster.names <- unique(Idents(SO.integrated))[order(unique(Idents(SO.integrated)))]
					markers.fm.list <- list()
						if(length(cluster.names)>1){ #do not calculate within group markers if group contains only one celltype
							for (i in 1:length(cluster.names)) {
									print(paste0("calculating markers for cluster ",cluster.names[i],". Total: ",length(cluster.names)," clusters"))
									markers.fm.list[[i]] <- FindMarkers(SO.integrated, ident.1 = cluster.names[i], min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
									markers.fm.list[[i]] <- markers.fm.list[[i]][which(markers.fm.list[[i]]$p_val_adj < (p.val.cutoff) ),] #select genes with p_val_adj > p.val.cutoff setting
							}	
						}
					IEmarkers[[j]] <- unlist(sapply(markers.fm.list, rownames))
					SO.integrated <- SO.backup
				}
				IEmarkers[[j+1]] <- "to be deleted" ;names(IEmarkers) <- c(names(clusters.type),"delete me");IEmarkers <- IEmarkers[1:j]   #necessary otherwise last element of list can be empty and will be deleted
			} 
			
		# Calculate music-pregrouped proportions
			Est.prop.bulk <- music_prop.cluster(bulk.eset = pseudobulk.eset, sc.eset = sc.eset, group.markers = IEmarkers, clusters = 'cellType', groups = 'clusterType', samples = 'SubjectName', clusters.type = clusters.type)
			estimate.music.pregrouped <- Est.prop.bulk$Est.prop.weighted.cluster

			correlation.result <- Eval_multi(prop.real = data.matrix(real.normalized), prop.est =  data.matrix(estimate.music.pregrouped), method.name = 'MuSiC pre-grouped')
			correlation.result

			#create comp fig heatmap as shown in MuSiC tutorial
				prop.comp.fig = Prop_comp_multi(prop.real = data.matrix(real.normalized),
												prop.est = 	data.matrix(estimate.music.pregrouped),
												method.name = c('MuSiC'), 
												title = 'Heatmap of Real and Est. Prop' )
}	

	#3. plot benchmark figures
		#creat correlation data set
			values   = c(as.numeric(real.normalized),as.numeric(estimate.music.pregrouped)) #values
			celltype = c(rep(colnames(real.normalized),each = nrow(real.normalized)),rep(colnames(estimate.music.pregrouped),each = nrow(estimate.music.pregrouped))) #celltype
			algorithm = c(rep("real.normalized", prod(dim(real.normalized))),rep("music", prod(dim(estimate.music.pregrouped)))) #algorithm
			correlation.data <- data.frame(values, celltype,algorithm)

		# correlation plot
			real.normalized.ordered <- real.normalized[,order(colnames(real.normalized))]
			estimate.music.pregrouped.ordered <- estimate.music.pregrouped[,order(colnames(estimate.music.pregrouped))]
			
			benchmark.df <- cbind(as.vector(real.normalized.ordered), as.vector(estimate.music.pregrouped.ordered),rep(colnames(real.normalized.ordered),each=5))
			colnames(benchmark.df) <- c("pseudobulk","MuSiC_estimate","celltype")
			benchmark.df <- as.data.frame(benchmark.df)
			benchmark.df[,1] <- as.numeric(benchmark.df[,1])
			benchmark.df[,2] <- as.numeric(benchmark.df[,2])
			benchmark.df[,3] <- as.factor(benchmark.df[,3])
				
			cor.val <- round(cor(benchmark.df$pseudobulk, benchmark.df$MuSiC_estimate), 2) #calculate pearson correlation
			cor.p <- cor.test(benchmark.df$pseudobulk, benchmark.df$MuSiC_estimate)$p.value
			
			#change colors
				color.vector <- pal_npg("nrc")(10) #requires ggsci package
				color.vector <- c(color.vector,"brown","khaki","orange3","orchid2")
			
			#change color order 
				benchmark.df$celltype <- factor(benchmark.df$celltype, levels = c(c("IgA plasma cell","IgG plasma cell","IgM plasma cell","B cell","cholangiocyte","stromal cell","hepatocyte","endothelial cell","T cell","NK CD56bright","NK CD56dim","monocyte","cDC","kupffer cell")))

			#visualize
			p=ggplot(benchmark.df, aes(x = pseudobulk, y = MuSiC_estimate, col=celltype)) + geom_point(alpha=0.6, size=3) + 
				geom_abline(color = "red", slope = 1) + 
				annotate(x = 0.05, y = 0.80, hjust=0,  geom = "text", label = paste0("Pearson's r","=",cor.val), size = 6) +
				annotate(x = 0.05, y = 0.75, hjust=0,  geom = "text", label = paste0("p-value=",signif(cor.p,3)), size = 6) +
				xlab("Real proportion (ground truth)") + ylab("MuSiC estimate") +
				theme_classic() +
				theme(axis.text.x = element_text(size = 20, color="black")) +
				theme(axis.text.y = element_text(size = 20, color="black")) +
				theme(axis.title.x = element_text(size = 20)) +
				theme(axis.title.y = element_text(size = 20)) +
				theme(legend.title = element_blank()) +
				theme(legend.text = element_text(size = 20))+
				scale_color_manual(values=color.vector)				
			p
	
		
# Part 3: perform pre-grouped MuSiC deconvolution on Training cohort 
#####################################################################

	#1. import training_56 data
		rawcount.file = "data_files/collapsed_data/GSE237330_training_raw_counts.gct"
		read.gct<-function(filename="NULL"){
			if (regexpr(".gct$",filename)==-1){
				stop("### input data should be .gct file! ###")
			}
			data<-read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T,check.names=F)
			data<-data[-1]
			return(data)
		}
			
		my.matrix      <- read.gct(file=rawcount.file)
		

	#2. create bulk eset
	{
		individual.ids <-  colnames(my.matrix)
		print(paste("there are", length(unique(individual.ids)),"unique individuals:", paste(unique(individual.ids),collapse=" , "))) #print as quality control
		
		base::names(individual.ids) <- colnames(my.matrix)
		individual.ids <- base::factor(individual.ids)
		n.individuals <- base::length(base::levels(individual.ids))
		sample.ids <- colnames(my.matrix)
		bulk.pheno <- base::data.frame(check.names=F, check.rows=F,stringsAsFactors=F,row.names=sample.ids,SubjectName=individual.ids,cellType=colnames(my.matrix))
		bulk.meta <- base::data.frame(labelDescription=base::c("SubjectName","cellType"),row.names=base::c("SubjectName","cellType"))
		bulk.pdata <- methods::new("AnnotatedDataFrame",data=bulk.pheno,varMetadata=bulk.meta)
		bulk.data <- base::as.matrix(my.matrix)
		bulk.eset <- Biobase::ExpressionSet(assayData=bulk.data,phenoData=bulk.pdata) # finale eset!! 
	}

	#3. calculate music-pregrouped proportions
		Est.prop.bulk <- music_prop.cluster(bulk.eset = bulk.eset, sc.eset = sc.eset, group.markers = IEmarkers, clusters = 'cellType', groups = 'clusterType', samples = 'SubjectName', clusters.type = clusters.type)
		estimate.music.pregrouped <- as.data.frame(Est.prop.bulk$Est.prop.weighted.cluster)
		my.celltypes <- colnames(estimate.music.pregrouped)

		
	#4. load NTP predicion
		{
		input.NTP.file = "training_NTP_prediction_result.txt"
		select.p = "BH.FDR"
		p.cutoff = 0.05
		
		#read NTP prediction file generated by NTP_KM.R 
			NTP.df <- read.table(file=input.NTP.file, header=T, sep="\t") 
			colnames(NTP.df)[1:2] <- c("SubjectName","PredictedCls") #Change first two colnames: "SubjectName","PredictedCls"
			rownames(NTP.df) <- NTP.df$SubjectName
			
			if(select.p == "BH.FDR"){NTP.df.intermediate    <- cbind(NTP.df$SubjectName,NTP.df$PredictedCls,NTP.df$BH.FDR, rep(0,nrow(NTP.df)))}
			if(select.p == "nominal.p"){NTP.df.intermediate <- cbind(NTP.df$SubjectName,NTP.df$PredictedCls,NTP.df$nominal.p, rep(0,nrow(NTP.df)))}
			for(i in 1:nrow(NTP.df.intermediate)){
				if(NTP.df.intermediate[i,3] < p.cutoff) { NTP.df.intermediate[i,4] <- NTP.df.intermediate[i,2] }else{ NTP.df.intermediate[i,4] <- "intermediate"}
				if(NTP.df.intermediate[i,4] == 1) { NTP.df.intermediate[i,4] <- "poor"}
				if(NTP.df.intermediate[i,4] == 2) { NTP.df.intermediate[i,4] <- "good"}
			}
			NTP.df$Prognosis <- NTP.df.intermediate[,4]
		}
		

		#read meta data & add prognosis
			training.meta.df <- read.table(file="data_files/meta_data/training_full_meta.txt",header=T)
			training.full.meta.prognosis <- merge(training.meta.df,NTP.df, by="SubjectName")
		
		#save full meta data with prognosis		
			#write.table(training.full.meta.prognosis, file="data_files/meta_data/training_full_meta_prognosis.txt", quote=F, row.names=F, sep="\t")  


	#5. create boxplots
		#create boxplot.df
			estimate.music.pregrouped <- as.matrix(estimate.music.pregrouped)
			training.meta.df <- training.full.meta.prognosis [, colnames(training.full.meta.prognosis ) %in% c("SubjectName","Prognosis","Metavir")]
			
			boxplot.df <- matrix(0,nrow(estimate.music.pregrouped)*ncol(estimate.music.pregrouped),3 + ncol(training.meta.df))
			colnames(boxplot.df) <- c("value","celltype","patient", colnames(training.meta.df))
			boxplot.df[,1] <- as.vector(estimate.music.pregrouped)
			boxplot.df[,2] <- rep(colnames(estimate.music.pregrouped), each = nrow(estimate.music.pregrouped))
			boxplot.df[,3] <- rep(rownames(estimate.music.pregrouped), ncol(estimate.music.pregrouped))
			
			for(i in 1:ncol(training.meta.df)){
				boxplot.df[,3+i] <- rep(training.meta.df[,i], ncol(estimate.music.pregrouped))
			}
			boxplot.df <- as.data.frame(boxplot.df)
			boxplot.df$value <- as.numeric(boxplot.df$value)

		#plot data
			#RISK-STRATIFIED boxplot prognosis Poor - intermediate - Good with individual datapoints, publication 
				ggplot(boxplot.df, aes(x=celltype, y=value, fill=Prognosis)) + 
					geom_boxplot(outlier.shape=NA) + 
					geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
					scale_fill_manual(values=c("lightblue","gray","tomato")) +
					xlab("") + ylab("estimated celltype proportion") +				
					theme_classic()+ 
					theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=.8, hjust=0.8)) 
			
				#set celltype comparisons for poor-intermediate-good 
					my_comparisons = list(c("poor", "good"))
					#my_comparisons = list(c("poor", "good"), c("poor", "intermediate"),c("good","intermediate"))
				
				#hepatocytes
					ph=ggplot(boxplot.df[boxplot.df$celltype=="hepatocyte",], aes(x=Prognosis, y=value, fill=Prognosis)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("lightblue","gray","tomato")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("hepatocyte")
					ph

				#endothelial cell
					pe=ggplot(boxplot.df[boxplot.df$celltype=="endothelial cell",], aes(x=Prognosis, y=value, fill=Prognosis)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("lightblue","gray","tomato")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("endothelial cell")
					pe
				
				#stromal cell
					ps=ggplot(boxplot.df[boxplot.df$celltype=="stromal cell",], aes(x=Prognosis, y=value, fill=Prognosis)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("lightblue","gray","tomato")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("stromal cell")
					ps
				
				#IgA plasma cell
					pa=ggplot(boxplot.df[boxplot.df$celltype=="IgA plasma cell",], aes(x=Prognosis, y=value, fill=Prognosis)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("lightblue","gray","tomato")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("IgA plasma cell")
					pa

				#combined
					library(patchwork)

					ph <- ph + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) 
					pe <- pe + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ylab("")
					ps <- ps + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ylab("")
					pa <- pa + ylab("") +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
					pc=wrap_plots(list(ph,pe,ps,pa), ncol=4)
					pc

			#METAVIR STRAIFIED box plot
				boxplot.df$Alt_metavir <- boxplot.df$Metavir 
				boxplot.df$Alt_metavir2 <- boxplot.df$Metavir 
				boxplot.df$Alt_metavir2[boxplot.df$Alt_metavir2 == "0_1" |boxplot.df$Alt_metavir2 == "2"] <- "0_1_2" #add column for 0_1 and 2 together

				ggplot(boxplot.df, aes(x=celltype, y=value, fill=Alt_metavir)) + 
					geom_boxplot(outlier.shape=NA) + 
					geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
					scale_fill_manual(values=c("peachpuff","tan1","tan3")) +
					xlab("") + ylab("estimated proportion") +				
					theme_classic()+ 
					theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=.8, hjust=0.8)) 
			
				#set celltype comparisons for Metavir
					my_comparisons = list(c("0_1", "3_4"))
					#my_comparisons = list(c("0_1", "3_4"), c("0_1", "2"),c("2", "3_4"))
					#my_comparisons = list(c("0_1", "3_4"),c("2", "3_4"))
					
					ph1=ggplot(boxplot.df[boxplot.df$celltype=="hepatocyte",], aes(x=Alt_metavir, y=value, fill=Alt_metavir)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("peachpuff","tan1","tan3"), labels=c("0 or 1","2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("hepatocyte")
					
					ph2=ggplot(boxplot.df[boxplot.df$celltype=="hepatocyte",], aes(x=Alt_metavir2, y=value, fill=Alt_metavir2)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("tan2","tan3"), labels=c("0 or 1 or 2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = list(c("0_1_2", "3_4")), size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("hepatocyte")
					ph=wrap_plots(list(ph1,ph2), ncol=1)

				#endothelial cell
					pe1=ggplot(boxplot.df[boxplot.df$celltype=="endothelial cell",], aes(x=Alt_metavir, y=value, fill=Alt_metavir)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("peachpuff","tan1","tan3"), labels=c("0 or 1","2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("endothelial cell")

					pe2=ggplot(boxplot.df[boxplot.df$celltype=="endothelial cell",], aes(x=Alt_metavir2, y=value, fill=Alt_metavir2)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("tan2","tan3"), labels=c("0 or 1 or 2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = list(c("0_1_2", "3_4")), size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("endothelial cell")
					pe=wrap_plots(list(pe1,pe2), ncol=1)
				
				#stromal cell
					ps1=ggplot(boxplot.df[boxplot.df$celltype=="stromal cell",], aes(x=Alt_metavir, y=value, fill=Alt_metavir)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("peachpuff","tan1","tan3"), labels=c("0 or 1","2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						guides(fill=guide_legend("Metavir")) +
						xlab("stromal cell")

					ps2=ggplot(boxplot.df[boxplot.df$celltype=="stromal cell",], aes(x=Alt_metavir2, y=value, fill=Alt_metavir2)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("tan2","tan3"), labels=c("0 or 1 or 2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = list(c("0_1_2", "3_4")), size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("stromal cell")
					ps=wrap_plots(list(ps1,ps2), ncol=1)
			
				#IgA plasma cell
					pa1=ggplot(boxplot.df[boxplot.df$celltype=="IgA plasma cell",], aes(x=Alt_metavir, y=value, fill=Alt_metavir)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("peachpuff","tan1","tan3"), labels=c("0 or 1","2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						guides(fill=guide_legend("Metavir")) +
						xlab("IgA plasma cell")

					pa2=ggplot(boxplot.df[boxplot.df$celltype=="IgA plasma cell",], aes(x=Alt_metavir2, y=value, fill=Alt_metavir2)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("tan2","tan3"), labels=c("0 or 1 or 2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = list(c("0_1_2", "3_4")), size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("IgA plasma cell") +
						guides(fill=guide_legend("Metavir")) 
					pa=wrap_plots(list(pa1,pa2), ncol=1)
											
				#combined
					library(patchwork)

					ph1 <- ph1 + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) 
					ph2 <- ph2 + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) 
					pe1 <- pe1 + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ylab("")
					pe2 <- pe2 + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ylab("")
					ps1 <- ps1 + ylab("") +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
					ps2 <- ps2 + ylab("") +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
					pc=wrap_plots(list(ph1,pe1,ps1,ph2,pe2,ps2), ncol=3)
					pc


# Part 4: perform pre-grouped MuSiC deconvolution on Validation cohort 
########################################################################
	#1. import validation data
		rawcount.file = "data_files/collapsed_data/GSE237331_validation_raw_counts.gct"
		read.gct<-function(filename="NULL"){
			if (regexpr(".gct$",filename)==-1){
				stop("### input data should be .gct file! ###")
			}
			data<-read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T,check.names=F)
			data<-data[-1]
			return(data)
		}
			
		my.matrix      <- read.gct(file=rawcount.file)
		

	#2. create bulk eset
	{
		individual.ids <-  colnames(my.matrix)
		print(paste("there are", length(unique(individual.ids)),"unique individuals:", paste(unique(individual.ids),collapse=" , "))) #print as quality control
		
		base::names(individual.ids) <- colnames(my.matrix)
		individual.ids <- base::factor(individual.ids)
		n.individuals <- base::length(base::levels(individual.ids))
		sample.ids <- colnames(my.matrix)
		bulk.pheno <- base::data.frame(check.names=F, check.rows=F,stringsAsFactors=F,row.names=sample.ids,SubjectName=individual.ids,cellType=colnames(my.matrix))
			
		bulk.meta <- base::data.frame(labelDescription=base::c("SubjectName","cellType"),row.names=base::c("SubjectName","cellType"))
		bulk.pdata <- methods::new("AnnotatedDataFrame",data=bulk.pheno,varMetadata=bulk.meta)
		bulk.data <- base::as.matrix(my.matrix)
		bulk.eset <- Biobase::ExpressionSet(assayData=bulk.data,phenoData=bulk.pdata) # finale eset!! 
	}

	#3. calculate music-pregrouped proportions
		Est.prop.bulk <- music_prop.cluster(bulk.eset = bulk.eset, sc.eset = sc.eset, group.markers = IEmarkers, clusters = 'cellType', groups = 'clusterType', samples = 'SubjectName', clusters.type = clusters.type)
		estimate.music.pregrouped <- as.data.frame(Est.prop.bulk$Est.prop.weighted.cluster)
		my.celltypes <- colnames(estimate.music.pregrouped)
		
	#4. load NTP predicion
		{
		input.NTP.file = "validation_NTP_prediction_result.txt"
		select.p = "BH.FDR"
		p.cutoff = 0.05
		
		#read NTP prediction file generated by NTP_KM.R 
			NTP.df <- read.table(file=input.NTP.file, header=T, sep="\t") 
			colnames(NTP.df)[1:2] <- c("SubjectName","PredictedCls") #Change first two colnames: "SubjectName","PredictedCls"
			rownames(NTP.df) <- NTP.df$SubjectName
			
			if(select.p == "BH.FDR"){NTP.df.intermediate    <- cbind(NTP.df$SubjectName,NTP.df$PredictedCls,NTP.df$BH.FDR, rep(0,nrow(NTP.df)))}
			if(select.p == "nominal.p"){NTP.df.intermediate <- cbind(NTP.df$SubjectName,NTP.df$PredictedCls,NTP.df$nominal.p, rep(0,nrow(NTP.df)))}
			for(i in 1:nrow(NTP.df.intermediate)){
				if(NTP.df.intermediate[i,3] < p.cutoff) { NTP.df.intermediate[i,4] <- NTP.df.intermediate[i,2] }else{ NTP.df.intermediate[i,4] <- "intermediate"}
				if(NTP.df.intermediate[i,4] == 1) { NTP.df.intermediate[i,4] <- "poor"}
				if(NTP.df.intermediate[i,4] == 2) { NTP.df.intermediate[i,4] <- "good"}
			}
			NTP.df$Prognosis <- NTP.df.intermediate[,4]
		}
		

		#read meta data & add prognosis
			validation.meta.df <- read.table(file="data_files/meta_data/validation_full_meta.txt",header=T)
			validation.full.meta.prognosis <- merge(validation.meta.df,NTP.df, by="SubjectName")
		
		#save full meta data with prognosis		
			#write.table(validation.full.meta.prognosis, file="data_files/meta_data/validation_full_meta_prognosis.txt", quote=F, row.names=F, sep="\t")  


	#5. create boxplots
		#create boxplot.df
			estimate.music.pregrouped <- as.matrix(estimate.music.pregrouped)
			validation.meta.df <- validation.full.meta.prognosis [, colnames(validation.full.meta.prognosis ) %in% c("SubjectName","Prognosis","Metavir")]
			
			boxplot.df <- matrix(0,nrow(estimate.music.pregrouped)*ncol(estimate.music.pregrouped),3 + ncol(validation.meta.df))
			colnames(boxplot.df) <- c("value","celltype","patient", colnames(validation.meta.df))
			boxplot.df[,1] <- as.vector(estimate.music.pregrouped)
			boxplot.df[,2] <- rep(colnames(estimate.music.pregrouped), each = nrow(estimate.music.pregrouped))
			boxplot.df[,3] <- rep(rownames(estimate.music.pregrouped), ncol(estimate.music.pregrouped))
			
			for(i in 1:ncol(validation.meta.df)){
				boxplot.df[,3+i] <- rep(validation.meta.df[,i], ncol(estimate.music.pregrouped))
			}
			
			boxplot.df <- as.data.frame(boxplot.df)
			boxplot.df$value <- as.numeric(boxplot.df$value)

		#plot data 
			#RISK-STRATIFIED boxplot prognosis Poor - intermediate - Good with individual datapoints, publication 
				ggplot(boxplot.df, aes(x=celltype, y=value, fill=Prognosis)) + 
					geom_boxplot(outlier.shape=NA) + 
					geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
					scale_fill_manual(values=c("lightblue","gray","tomato")) +
					xlab("") + ylab("estimated celltype proportion") +				
					theme_classic()+ 
					theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=.8, hjust=0.8)) 
			
				#set celltype comparisons for poor-intermediate-good 
					my_comparisons = list(c("poor", "good"))
					#my_comparisons = list(c("poor", "good"), c("poor", "intermediate"),c("good","intermediate"))
				
				#hepatocytes
					ph=ggplot(boxplot.df[boxplot.df$celltype=="hepatocyte",], aes(x=Prognosis, y=value, fill=Prognosis)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("lightblue","gray","tomato")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("hepatocyte")
					ph

				#endothelial cell
					pe=ggplot(boxplot.df[boxplot.df$celltype=="endothelial cell",], aes(x=Prognosis, y=value, fill=Prognosis)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("lightblue","gray","tomato")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("endothelial cell")
					pe
				
				#stromal cell
					ps=ggplot(boxplot.df[boxplot.df$celltype=="stromal cell",], aes(x=Prognosis, y=value, fill=Prognosis)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("lightblue","gray","tomato")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("stromal cell")
					ps
				
				#IgA plasma cell
					pa=ggplot(boxplot.df[boxplot.df$celltype=="IgA plasma cell",], aes(x=Prognosis, y=value, fill=Prognosis)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("lightblue","gray","tomato")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("IgA plasma cell")
					pa

				#combined
					library(patchwork)

					ph <- ph + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) 
					pe <- pe + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ylab("")
					ps <- ps + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ylab("")
					pa <- pa + ylab("") +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
					pc=wrap_plots(list(ph,pe,ps,pa), ncol=4)
					pc
					
			#METAVIR STRATIFIED box plot
				boxplot.df$Alt_metavir <- boxplot.df$Metavir 
				boxplot.df$Alt_metavir2 <- boxplot.df$Metavir 
				boxplot.df$Alt_metavir2[boxplot.df$Alt_metavir2 == "0_1" |boxplot.df$Alt_metavir2 == "2"] <- "0_1_2" #add column for 0_1 and 2 together

				ggplot(boxplot.df, aes(x=celltype, y=value, fill=Alt_metavir)) + 
					geom_boxplot(outlier.shape=NA) + 
					geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
					scale_fill_manual(values=c("peachpuff","tan1","tan3")) +
					xlab("") + ylab("estimated proportion") +				
					theme_classic()+ 
					theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=.8, hjust=0.8)) 
			
				#set celltype comparisons for Metavir
					my_comparisons = list(c("0_1", "3_4"))
					#my_comparisons = list(c("0_1", "3_4"), c("0_1", "2"),c("2", "3_4"))
					#my_comparisons = list(c("0_1", "3_4"),c("2", "3_4"))
					
					ph1=ggplot(boxplot.df[boxplot.df$celltype=="hepatocyte",], aes(x=Alt_metavir, y=value, fill=Alt_metavir)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("peachpuff","tan1","tan3"), labels=c("0 or 1","2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("hepatocyte")
					
					ph2=ggplot(boxplot.df[boxplot.df$celltype=="hepatocyte",], aes(x=Alt_metavir2, y=value, fill=Alt_metavir2)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("tan2","tan3"), labels=c("0 or 1 or 2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = list(c("0_1_2", "3_4")), size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("hepatocyte")
					ph=wrap_plots(list(ph1,ph2), ncol=1)

				#endothelial cell
					pe1=ggplot(boxplot.df[boxplot.df$celltype=="endothelial cell",], aes(x=Alt_metavir, y=value, fill=Alt_metavir)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("peachpuff","tan1","tan3"), labels=c("0 or 1","2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("endothelial cell")

					pe2=ggplot(boxplot.df[boxplot.df$celltype=="endothelial cell",], aes(x=Alt_metavir2, y=value, fill=Alt_metavir2)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("tan2","tan3"), labels=c("0 or 1 or 2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = list(c("0_1_2", "3_4")), size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("endothelial cell")
					pe=wrap_plots(list(pe1,pe2), ncol=1)
				
				#stromal cell
					ps1=ggplot(boxplot.df[boxplot.df$celltype=="stromal cell",], aes(x=Alt_metavir, y=value, fill=Alt_metavir)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("peachpuff","tan1","tan3"), labels=c("0 or 1","2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						guides(fill=guide_legend("Metavir")) +
						xlab("stromal cell")

					ps2=ggplot(boxplot.df[boxplot.df$celltype=="stromal cell",], aes(x=Alt_metavir2, y=value, fill=Alt_metavir2)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("tan2","tan3"), labels=c("0 or 1 or 2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = list(c("0_1_2", "3_4")), size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("stromal cell")
					ps=wrap_plots(list(ps1,ps2), ncol=1)
			
				#IgA plasma cell
					pa1=ggplot(boxplot.df[boxplot.df$celltype=="IgA plasma cell",], aes(x=Alt_metavir, y=value, fill=Alt_metavir)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("peachpuff","tan1","tan3"), labels=c("0 or 1","2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						guides(fill=guide_legend("Metavir")) +
						xlab("IgA plasma cell")

					pa2=ggplot(boxplot.df[boxplot.df$celltype=="IgA plasma cell",], aes(x=Alt_metavir2, y=value, fill=Alt_metavir2)) + 
						geom_boxplot(outlier.shape=NA) + 
						geom_point(aes(color=Prognosis,size = 2), shape = 21,size = 2,colour = "black", position=position_jitterdodge(jitter.width=0.05)) +
						scale_fill_manual(values=c("tan2","tan3"), labels=c("0 or 1 or 2","3 or 4")) +
						ylab("estimated proportion") +	
						theme_classic()+ 
						stat_compare_means(method = "wilcox.test", comparisons = list(c("0_1_2", "3_4")), size=5) +
						theme(axis.text.x = element_blank()) + 
						theme(axis.title.x = element_text(size=15, color="black")) +
						theme(axis.text.y = element_text(size=15, color="black")) +
						theme(axis.ticks.x = element_blank()) +
						theme(axis.title.y = element_text(size=15, color="black")) +
						theme(legend.title = element_text(size=15, color="black")) +
						theme(legend.text = element_text(size=15, color="black")) +
						xlab("IgA plasma cell") +
						guides(fill=guide_legend("Metavir")) 
					pa=wrap_plots(list(pa1,pa2), ncol=1)
											
				#combined
					library(patchwork)

					ph1 <- ph1 + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) 
					ph2 <- ph2 + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) 
					pe1 <- pe1 + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ylab("")
					pe2 <- pe2 + theme(legend.position='none') +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ylab("")
					ps1 <- ps1 + ylab("") +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
					ps2 <- ps2 + ylab("") +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
					pc=wrap_plots(list(ph1,pe1,ps1,ph2,pe2,ps2), ncol=3)
					pc
	