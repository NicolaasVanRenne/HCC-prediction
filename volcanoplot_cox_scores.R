##############################
# Volcano plot of cox scores
# of 557 gene signature
##############################
#load library
	library(ggplot2)
	library(ggrepel)
	library(dplyr)

#set working directory
	setwd("C:/my_dir")
	
#load data
	volcano.df <- read.delim(file="training_RPM_filtered_signature.txt",header=T,check.names=F)		
	volcano.df <- as.data.frame(volcano.df)
	rownames(volcano.df) <- volcano.df$probeid
	
#Genes positively or negatively correlated with HCC-incidence were selected using the cox score
	ggplot(data=volcano.df, aes(x=statistic, y=-log10(nominal.p))) + geom_point() + theme_classic() 
	
#Publication quality volcano plot
	#load poor and good prognosis genes
		signature.df <- read.delim(file="training_RPM_filtered_signature_result_ordered.txt",header=T,check.names=F)		
		poor.genes <- signature.df[signature.df$statistic > 0, 1]
		good.genes <- signature.df[signature.df$statistic < 0, 1]
		
	#add column of colored dots based on risk genes 
		volcano.df$Direction <- "None"
		volcano.df$Direction[rownames(volcano.df) %in% poor.genes] <- "HCC-incidence correlated" 		
		volcano.df$Direction[rownames(volcano.df) %in% good.genes] <- "HCC-incidence anti-correlated"			
		
		volcano.df$Direction <- as.factor(volcano.df$Direction)
		volcano.df$Direction <- factor(volcano.df$Direction, levels = c("HCC-incidence anti-correlated", "None", "HCC-incidence correlated"))
		ggplot(data=volcano.df, aes(x=statistic, y=-log10(nominal.p), col=Direction)) + geom_point() +  theme_classic() + 	scale_color_manual(values=c("blue", "black", "red")) 
	
		# partially transparent points by setting `alpha`
			linecolors <- c("navy", "darkgrey", "firebrick")
			fillcolors <- c("navy", "darkgrey", "firebrick")
			
		#set plot order
			plot.order <- rep(0,nrow(volcano.df))
			plot.order[volcano.df$Direction =="None"] <- 1
			plot.order[volcano.df$Direction =="HCC-incidence anti-correlated"] <- 2
			plot.order[volcano.df$Direction =="HCC-incidence correlated"] <- 3
			volcano.df$plotorder <- plot.order
		
		#log10 pval
		volcano.log10.df <- volcano.df
		volcano.log10.df$nominal.p[volcano.log10.df$nominal.p == 0] <- 0.001 #change pval 0 to 0.001 
		volcano.log10.df$log10pval <- -log(volcano.log10.df$nominal.p,10)
		
		ggplot(volcano.log10.df %>% arrange(plotorder),  aes(x=statistic, y=log10pval, col=Direction, fill=Direction)) + 
			geom_point(shape = 21, alpha = 0.2, size = 1) +
			scale_color_manual(values=linecolors) +
			scale_fill_manual(values=fillcolors) +
			geom_hline(yintercept=-0.1, linetype='dashed', col = 'darkgrey') +
			labs(y= "P value", x = "cox score") +
			theme_classic()


	#add labels
		library(ggrepel)

		#select only IGHA1 and IGHA2 label
			options(ggrepel.max.overlaps = Inf)
			UP.genes <- c("IGHA2","IGHA1")
			DN.genes <- c() #c() 
					
			volcano.log10.df$genelabel <- NA
			volcano.log10.df$genelabel[rownames(volcano.log10.df) %in% UP.genes] <- rownames(volcano.log10.df)[rownames(volcano.log10.df) %in% UP.genes]
			volcano.log10.df$genelabel[rownames(volcano.log10.df) %in% DN.genes] <- rownames(volcano.log10.df)[rownames(volcano.log10.df) %in% DN.genes]

			p=ggplot(volcano.log10.df %>% arrange(plotorder),  aes(x=statistic, y=log10pval, col=Direction, fill=Direction, label=genelabel)) + 
				geom_point(shape = 21, alpha = 0.2, size = 1) +
				xlim(c(-3,8))+
				scale_color_manual(values=linecolors) +
				scale_fill_manual(values=fillcolors) +
				geom_hline(yintercept=1, linetype='dashed', col = 'darkgrey') +
				labs(y= "-log10(p-value)", x = "cox score") +
				theme_classic()	+
				theme(axis.text.x = element_text(size = 20, color="black")) +
				theme(axis.title.x = element_text(size = 20)) +
				theme(axis.text.y = element_text(size = 20, color="black")) +
				theme(axis.title.y = element_text(size = 20)) +
				theme(legend.position="top", legend.text = element_text(size=10), legend.title= element_text(size=10)) +guides(col=guide_legend(override.aes = list(size=5), nrow=3)) +
				geom_text_repel(
					size   = 5,
					force        = 1,
					nudge_x      = 10,
					direction    = "y",
					segment.size = 0.2, #	line segment thickness
					show.legend=FALSE
			) 
			p
				
										