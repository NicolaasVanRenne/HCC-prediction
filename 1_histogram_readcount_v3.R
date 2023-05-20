############################################
# 1. create figures of readcount histograms
############################################	
	library(ggplot2)
	library(patchwork)
	
	save.plots=FALSE
	
	#set working directory
		setwd("C:/my_directory")

	#create read.gct function
		read.gct<-function(filename){
			data<-read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T,check.names=F)
			data<-data[-1]
			return(data)
		}

	#read training raw readcount file
		training.data <- read.gct("data_files/processed_data/training_56_raw.gct")
		
	#read validation raw readcount file
		validation.data <- read.gct("data_files/processed_data/validation_51_raw.gct")
	
	
	#calculate readcounts
		training.reads <- data.frame(colSums(training.data)/10^6)
		colnames(training.reads) <- "readcount"
		
		validation.reads <- data.frame(colSums(validation.data)/10^6)
		colnames(validation.reads) <- "readcount"		
		
	#create plots
		p1=ggplot(training.reads, aes(x=readcount)) + 
			geom_histogram(binwidth=1, center=0.5, colour="black", fill="grey") + 
			theme_classic() +
			scale_y_continuous(expand = c(0, 0), limits=c(0,8.5), breaks=c(0,2,4,6,8)) +
			labs(y="samples (n)", x="read count (million)") + 
			#ggtitle("Training set") + theme(plot.title = element_text(hjust=0.5, size = 15, face = "bold"))+
			theme(axis.text.x = element_text(size = 15, color="black"), axis.title.x = element_text(size = 15, color="black")) +
			theme(axis.text.y = element_text(size = 15, color="black"), axis.title.y = element_text(size = 15, color="black", margin=margin(r=10))) 
		p1

		p2=ggplot(validation.reads, aes(x=readcount)) + 
			geom_histogram(binwidth=1, center=0.5, colour="black", fill="grey") + 
			theme_classic() +
			scale_y_continuous(expand = c(0, 0), limits=c(0,6.5), breaks=c(0,2,4,6)) +
			scale_x_continuous(expand = c(0, 0), limits=c(-0.5,18), breaks=c(0,5,10,15)) +
			labs(y="samples (n)", x="read count (million)") + 
			#ggtitle("Training set") + theme(plot.title = element_text(hjust=0.5, size = 15, face = "bold"))+
			theme(axis.text.x = element_text(size = 15, color="black"), axis.title.x = element_text(size = 15, color="black")) +
			theme(axis.text.y = element_text(size = 15, color="black"), axis.title.y = element_text(size = 15, color="black", margin=margin(r=10))) 
		p2	
		
		p1+p2
	
	#save plots
		if(save.plots==TRUE){ggsave(p1, file= "output/plots/raw_readcount_training.tiff", width = 4.5, height=3)}
		if(save.plots==TRUE){ggsave(p2, file= "output/plots/raw_readcount_validation.tiff", width = 3, height=3)}
							
							
							
							