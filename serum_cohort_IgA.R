###########################################################
#
# analysis Serum cohort
#
##########################################################

save.plots=FALSE

library(ggplot2)
library(dplyr)
library(lubridate)
library(patchwork)
library(survival)
library(survminer)

#set dir
	setwd("C:/my_dir")

#load databases
	IgA.df 	<- read.table(file="data_files/serum_data/serum_cohort.txt", sep="\t", header=T) 

#preprocess data
	{
	#add fibrosis bin F2 = 0; F3 = 1; F4 = 2
		fibrosis.bin <- rep(0, nrow(IgA.df))
		for(i in seq_along(fibrosis.bin)){
			if(IgA.df$Fibrosis[i] == "F2")	 { fibrosis.bin[i] <-0 }
			if(IgA.df$Fibrosis[i] == "F3")	 { fibrosis.bin[i] <-1 }
			if(IgA.df$Fibrosis[i] == "F4")	 { fibrosis.bin[i] <-2 }
		}
		IgA.df$fibrosis.bin <- fibrosis.bin

	#transform censored time in years
		transform.censoredtime= "yes" #"yes" or "no"
		if(transform.censoredtime == "yes"){
			IgA.df$CensoredTime <- (IgA.df$CensoredTime / 365.24)
		}
	
	#add HBV
		IgA.df$HBV <- rep(NA,nrow(IgA.df))
		for(i in 1:nrow(IgA.df)){
			if(IgA.df$Virus[i] == "HBV"){ IgA.df$HBV[i] <- "yes" }
			if(IgA.df$Virus[i] != "HBV"){ IgA.df$HBV[i] <- "no" }
		}
		
	#add HCV
		IgA.df$HCV <- rep(NA,nrow(IgA.df))
		for(i in 1:nrow(IgA.df)){
			if(IgA.df$Virus[i] == "HCV"){ IgA.df$HCV[i] <- "yes" }
			if(IgA.df$Virus[i] != "HCV"){ IgA.df$HCV[i] <- "no" }
		}
		
	#add HBVHVC
		IgA.df$HBVHCV <- rep(NA,nrow(IgA.df))
		for(i in 1:nrow(IgA.df)){
			if(IgA.df$Virus[i] == "HBV HCV"){ IgA.df$HBVHCV[i] <- "yes" }
			if(IgA.df$Virus[i] != "HBV HCV"){ IgA.df$HBVHCV[i] <- "no" }
		}
		
		
	#add HBVDelta
		IgA.df$HBVDelta <- rep(NA,nrow(IgA.df))
		for(i in 1:nrow(IgA.df)){
			if(IgA.df$Virus[i] == "HBV-HDV"){ IgA.df$HBVDelta[i] <- "yes" }
			if(IgA.df$Virus[i] != "HBV-HDV"){ IgA.df$HBVDelta[i] <- "no" }
		}
	}

#patient table
	#number of patients:
		nrow(IgA.df)
	#Age median, IQR 
		quantile(IgA.df$Age)
	#serum IgA, IQR 
		quantile(IgA.df$SerumIgA)
	#follow-up, IQR
		quantile(IgA.df$CensoredTime)
	#HCC number
		sum(IgA.df$BinaryHCC)
		100*table(IgA.df$BinaryHCC)/nrow(IgA.df)
	#Virus
		#HCV
			table(IgA.df$HCV)
			100*table(IgA.df$HCV)/nrow(IgA.df)
		#HBV
			table(IgA.df$HBV)
			100*table(IgA.df$HBV)/nrow(IgA.df)
		#HBV-HCV
			table(IgA.df$HBVHCV)
			100*table(IgA.df$HBVHCV)/nrow(IgA.df)
		#HBV-Delta
			table(IgA.df$HBVDelta)
			100*table(IgA.df$HBVDelta)/nrow(IgA.df)

		#Gender
			table(IgA.df$Gender)
			100*table(IgA.df$Gender)/nrow(IgA.df)
		#Fibrosis
			table(IgA.df$fibrosis.bin)
			100*table(IgA.df$fibrosis.bin)/nrow(IgA.df)
		#Diagnostic
			table(IgA.df$Diagnostic)
			100*table(IgA.df$Diagnostic)/nrow(IgA.df)



# Cox regression analysis 
{		
	#Univariate cox regression
		#Age 
			coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ Age, data = IgA.df, method = "efron", robust = F)
			summary(coxph.result)
			wald.score.age 	<- as.numeric(coxph.result$wald.test) #wald test score
			wald.p.age 		<- as.numeric(summary(coxph.result)$waldtest[3]) #wald test score p
			
		#Serum IgA
			coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ SerumIgA, data = IgA.df, method = "efron", robust = F)
			summary(coxph.result)
			wald.score.serumIgA 	<- as.numeric(coxph.result$wald.test) #wald test score
			wald.p.serumIgA 		<- as.numeric(summary(coxph.result)$waldtest[3]) #wald test score p

		#Fibrosis coxph
			coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ fibrosis.bin, data = IgA.df, method = "efron", robust = F)
			summary(coxph.result)
			wald.score.fibrosis 	<- as.numeric(coxph.result$wald.test) #wald test score
			wald.p.fibrosis 		<- as.numeric(summary(coxph.result)$waldtest[3]) #wald test score p
		
		#HIV
			coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ HIV, data = IgA.df, method = "efron", robust = F)
			summary(coxph.result)
			wald.score.HIV 	<- as.numeric(coxph.result$wald.test) #wald test score
			wald.p.HIV 		<- as.numeric(summary(coxph.result)$waldtest[3]) #wald test score p

		#HCV
			coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ HCV, data = IgA.df, method = "efron", robust = F)
			summary(coxph.result)
			wald.score.HCV 	<- as.numeric(coxph.result$wald.test) #wald test score
			wald.p.HCV 		<- as.numeric(summary(coxph.result)$waldtest[3]) #wald test score p
	
		#HBV
			coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ HBV, data = IgA.df, method = "efron", robust = F)
			summary(coxph.result)
			wald.score.HBV 	<- as.numeric(coxph.result$wald.test) #wald test score
			wald.p.HBV 		<- as.numeric(summary(coxph.result)$waldtest[3]) #wald test score p
	
		#HBV-HCV
			coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ HBVHCV, data = IgA.df, method = "efron", robust = F)
			summary(coxph.result)
			wald.score.HBVHCV 	<- as.numeric(coxph.result$wald.test) #wald test score
			wald.p.HBVHCV 		<- as.numeric(summary(coxph.result)$waldtest[3]) #wald test score p
	
		#HBV-delta
			coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ HBVDelta, data = IgA.df, method = "efron", robust = F)
			summary(coxph.result)
			wald.score.HBVDelta 	<- as.numeric(coxph.result$wald.test) #wald test score
			wald.p.HBVDelta 		<- as.numeric(summary(coxph.result)$waldtest[3]) #wald test score p
	
		#Gender
			coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ Gender, data = IgA.df, method = "efron", robust = F)
			summary(coxph.result)
			wald.score.gender 	<- as.numeric(coxph.result$wald.test) #wald test score
			wald.p.gender 		<- as.numeric(summary(coxph.result)$waldtest[3]) #wald test score p
			
			#Create bar graph for Univariate Wald test scores
				library(ggplot2)
				data <- data.frame(
				Variable = c("Gender", "HBV", "HCV", "HBV-HDV", "HIV", "Age", "Fibrosis",  "Serum IgA"),
				Wald_Test_Score = c(wald.score.gender, wald.score.HBV, wald.score.HCV, wald.score.HBVDelta, wald.score.HIV, wald.score.age, wald.score.fibrosis,  wald.score.serumIgA),
				Wald_Test_P = c(wald.p.gender, wald.p.HBV, wald.p.HCV, wald.p.HBVDelta, wald.p.HIV, wald.p.age, wald.p.fibrosis,  wald.p.serumIgA)
				)
				data$P_value_label <- ifelse(data$Wald_Test_P >= 0.05, "NS", round(data$Wald_Test_P,6))
				
			# Create the bar chart HORIZONTAL
				data$Variable <- factor(data$Variable, levels = rev(data$Variable))
				p=ggplot(data, aes(x = Wald_Test_Score, y = Variable)) +
					geom_bar(stat = "identity", fill = "lightgrey", color = "black", width=0.8) +
					#geom_text(aes(label = P_value_label, x = Wald_Test_Score + 0.2),hjust = 0, size = 8, color = "black")+
					geom_text(aes(label = P_value_label, x = Wald_Test_Score + 2))+
					labs(title = NULL, x = "Wald Test Score (univariate)", y = NULL) +
					theme_classic() +
					theme(
						axis.title.x = element_text(size = 12, color="black"),
						axis.text.x = element_text(size = 12, color="black"),
						axis.text.y = element_text(size = 12, color="black"),
						plot.title = element_text(size = 12, color="black", hjust = 0.5)
					) +	
					scale_x_continuous(limits = c(0, max(data$Wald_Test_Score)+5), expand = c(0, 0))  # Ensure bars start at the origin
					#scale_x_continuous(expand = c(0, 0))  # Ensure bars start at the origin
				p
					if(save.plots==TRUE){  ggsave(p, file= "output/plots/wald_score_univariate_F2F3F4.tiff", width = 5, height=2)	}	
					if(save.plots==TRUE){  ggsave(p, file= "output/plots/wald_score_univariate_F2F3F4_large.tiff", width = 4, height=4)	}	
					
	#Multivariate cox regression
		#SerumIgA + Age + Fibrosis coxph
			coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ SerumIgA + Age + fibrosis.bin, data = IgA.df, method = "efron", robust = F)
			summary(coxph.result)
		
		#SerumIgA + Fibrosis coxph
			coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ SerumIgA + fibrosis.bin, data = IgA.df, method = "efron", robust = F)
			summary(coxph.result)
			summary(coxph.result)$waldtest[3]
}


#Akaike Information Criterion 
	coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ fibrosis.bin, data = IgA.df, method = "efron", robust = F)
	coxph.result ;AIC(coxph.result)
	coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ fibrosis.bin + SerumIgA  , data = IgA.df, method = "efron", robust = F)
	coxph.result ;AIC(coxph.result)
	
	
#AUROC
	library(pROC)
	
	#Single AUROCs
		#1. Serum IgA 
			roc.IgA <- roc(IgA.df$BinaryHCC, IgA.df$SerumIgA)
			auc(roc.IgA)
			plot(roc.IgA)		
		
		#2. Fibrosis.bin
			roc.fibrosis <- roc(IgA.df$BinaryHCC, IgA.df$fibrosis.bin)
			auc(roc.fibrosis)
			plot(roc.fibrosis)			
		
	
	#Combined AUROCs
		#Serum IgA + Fibrosis.bin
			coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ SerumIgA + fibrosis.bin  , data = IgA.df, method = "efron", robust = F)
			prediction <- predict(coxph.result)
			category <- IgA.df$BinaryHCC
			roc.sf <- roc(category, prediction)
			auc(roc.sf)
			plot(roc.sf)		


		# Plotting both ROC curves on the same plot
			plot(roc.IgA, col = "black", main = "ROC for Serum IgA and fibrosis")
			plot(roc.fibrosis, col = "purple", add = TRUE)
			plot(roc.sf, col = "orange", add = TRUE)
			# Adding legend
			legend("bottomright", legend = c("Serum IgA", "Fibrosis", "Serum IgA + Fibrosis"), col = c("black", "purple", "orange"), lwd = 2)

		#manual ROC cutoffs for IgA
			ROC.df <- IgA.df[order(IgA.df$SerumIgA),][,colnames(IgA.df) %in% c("SerumIgA","BinaryHCC")]
			n.events <- sum(ROC.df$BinaryHCC)
	
			for (i in 1:nrow(ROC.df)){
				ROC.df$IgAHCC[i] <- sum(ROC.df[ROC.df$SerumIgA == ROC.df$SerumIgA[i],]$BinaryHCC)
			}
			
			ROC.df <- ROC.df[!duplicated(as.numeric(ROC.df$SerumIgA)),] #remove duplicate values
			roc.result <- roc(IgA.df$BinaryHCC, IgA.df$SerumIgA) #roc(category,predictor)
			
			AUROC.table <- cbind(rbind(rep(0,2),ROC.df),roc.result$sensitivities,roc.result$specificities)
			AUROC.table
			
			
#Create Kaplan Meier curves Serum IgA only
	#set new cutoff serum IgA, cutoffs based on AUROC curves 
		cutofflow = 150 #serumIgA cutoff, lower or equal than this is low SerumIgACls
		cutoffhigh = 400 #serumIgA cutoff, higher or equal than this is high SerumIgACls 
		
	#add stratification to data 
		status <- rep(0, nrow(IgA.df))
		SerumIgACls <- rep(NA, nrow(IgA.df))
			for(i in 1:nrow(IgA.df)){ if(IgA.df$SerumIgA[i] >= cutoffhigh){ SerumIgACls[i] <- "high" }}
			for(i in 1:nrow(IgA.df)){ if(IgA.df$SerumIgA[i] <= cutofflow){ SerumIgACls[i] <- "low" }}
			for(i in 1:nrow(IgA.df)){ if(IgA.df$SerumIgA[i] < cutoffhigh & IgA.df$SerumIgA[i] > cutofflow){ SerumIgACls[i] <- "intermediate" }}
		IgA.df$SerumIgACls <- SerumIgACls

	#read clinical data		  
		km.df <- cbind(IgA.df$BinaryHCC, IgA.df$CensoredTime, IgA.df$SerumIgACls)
		colnames(km.df) <- c("HCC","time","SerumIgA")
		km.df <- as.data.frame(km.df)
	
		set.time  = "time"    #column name of time
		set.event = "HCC"   #column name 
		
	#run survfit
		survfit.result       				<- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ SerumIgA, data = km.df)
		survfit.intermediate.removed.result	<- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ SerumIgA, data = km.df[which(km.df$SerumIgA != "intermediate"),])
		survfit.low.removed.result			<- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ SerumIgA, data = km.df[which(km.df$SerumIgA != "low"),])
		survfit.high.removed.result			<- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ SerumIgA, data = km.df[which(km.df$SerumIgA != "high"),])

	#create truncated Kaplan Meier curve when patients at risk drops below 10% for visualization purpose (not for p-value!)
		#store p-values
			p.pairwise 			<- pairwise_survdiff(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ SerumIgA, data = km.df, p.adjust.method = "none")
			pval.low.vs.high 	<- signif(surv_pvalue(survfit.intermediate.removed.result)$pval ,3)
			pval.int.vs.high 	<- signif(surv_pvalue(survfit.low.removed.result)$pval ,3)
			pval.low.vs.int 	<- signif(surv_pvalue(survfit.high.removed.result)$pval ,3)

		#settings
			cutoff.at.risk.percentage = 10 #set percentage to cut eg. 'cutoff.at.risk = 10' for 10%
			#set.event = "HCC_event" #set event to find 10% at risk
		
		#calculate time for cutoff 
			km.truncate.df <-km.df
			km.truncate.df$all <- rep("all", nrow(km.df))
			km.truncate.df[,colnames(km.truncate.df) %in% set.event] <- rep(1,nrow(km.truncate.df))
			
			survfit.truncate <- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ all, data = km.truncate.df)
			summary(survfit.truncate) #show at time and risk table
		
			n.at.risk <- max(summary(survfit.truncate)$n.risk)
			n.at.risk.cutoff <- ceiling((n.at.risk/100)*cutoff.at.risk.percentage)
			cutoff.time.x.axis <- min(summary(survfit.truncate)$time[which(summary(survfit.truncate)$n.risk <= n.at.risk.cutoff)])
				{
					print(paste(n.at.risk,"patients in kaplan meier curve"))
					print(paste("This algorithm will cut off the kaplan meier curve when less than",cutoff.at.risk.percentage,"% of patients are at risk."))
					print(paste("This KM curve will be cut when less than",n.at.risk.cutoff,"patients are at risk."))
					print(paste("This KM curve will be cut at timepoint",cutoff.time.x.axis))
				}
	
				p=ggsurvplot(survfit.result,
					fun= "event",
					pval = paste0("high VS low p=",pval.low.vs.high,"\n", "high VS int. p=",pval.int.vs.high,"\n" , "int. VS low p=",pval.low.vs.int ),   
					pval.coord = c(0,0.50), 
					pval.size=7,
					conf.int = FALSE,
					risk.table = TRUE, risk.table.fontsize=6,
					risk.table.height=0.3,
					risk.table.y.text = TRUE,
					font.x=20,
					font.y=20,
					font.tickslab=20,
					palette = c("tomato", "darkgrey","mediumblue"),
					xlim = c(0,cutoff.time.x.axis), 
					break.time.by=3,
					legend = 'none',
					xlab="Time (years)",
					ylab="HCC incidence"
					)
					p$table <- p$table + 
							theme(axis.text.x = element_text(size = 20)) +
							theme(axis.title.x = element_text(size = 20)) +
							theme(axis.title.y = element_text(size = 20)) + 
							ylab("HCC risk") + 
							scale_y_discrete(labels=c('low', 'int.', 'high'))
					p$plot <- p$plot + geom_vline(xintercept = cutoff.time.x.axis, linetype = "dashed", color = "red", size = 1)# Modify the main plot to include a vertical line
				p
				if(save.plots==TRUE){  tiff('output/plots/KM_serumIgA_n106.tiff', units="in", width=5, height=7, res=300, compression = 'lzw');print(p);dev.off()  }
	
				p=ggsurvplot(survfit.result,
					fun= "event",
					#pval = paste0("high VS low p=",pval.low.vs.high,"\n", "high VS int. p=",pval.int.vs.high,"\n" , "int. VS low p=",pval.low.vs.int ),   
					#pval.coord = c(0,0.90), 
					#pval.size=7,
					conf.int = FALSE,
					risk.table = TRUE, risk.table.fontsize=6,
					risk.table.height=0.3,
					risk.table.y.text = TRUE,
					font.x=20,
					font.y=20,
					font.tickslab=20,
					palette = c("tomato", "darkgrey","mediumblue"),
					xlim = c(0,cutoff.time.x.axis), 
					break.time.by=3,
					legend = 'none',
					xlab="Time (years)",
					ylab="HCC incidence"
					)
					p$table <- p$table + 
							theme(axis.text.x = element_text(size = 20)) +
							theme(axis.title.x = element_text(size = 20)) +
							theme(axis.title.y = element_text(size = 20)) + 
							ylab("HCC risk") + 
							scale_y_discrete(labels=c('low', 'int.', 'high'))
				p
				if(save.plots==TRUE){  tiff('output/plots/KM_serumIgA_noP_n106.tiff', units="in", width=5, height=7, res=300, compression = 'lzw');print(p);dev.off()  }
	
	#FIND HCC INCIDENCE AT 3 YEARS
		summary(survfit.result, 3)
	



