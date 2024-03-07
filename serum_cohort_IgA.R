###########################################################
#
# analysis Serum cohort
#
##########################################################

#set working directory
	setwd("C:/my_dir")

#settings
	remove.1year = FALSE #TRUE or FALSE. TRUE if you want to removed all patients with less than 1 year censored time from the serum cohort

#load libraries
	library(ggplot2)
	library(dplyr)
	library(lubridate)
	library(patchwork)
	library(survival)
	library(survminer)

#load databases
	IgA.df 	<- read.table(file="data_files/serum_data/serum_cohort_binned.txt", sep="\t", header=T)

	if(remove.1year==TRUE){
		(print("removed all patients with less than 1 year censored time from the serum cohort, n= 169"))
		censoredtime.cutoff = 365.24 #set to zero to ignore
		IgA.df <- IgA.df[!IgA.df$CensoredTime < censoredtime.cutoff, ]
		nrow(IgA.df)
	}else{print("all patients of the serum cohort are used (n=188)")}
	
#add serum IgA binned per 50 mg/dL (0-49 = 0, 50-99 = 1, 100-149 = 2 etc..)
	IgA.df$SerumIgA50 <- floor(IgA.df$SerumIgA / 50)

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

#AUROC serum IgA (for n=188)
	if(remove.1year==TRUE){print("set remove.1year = FALSE to reproduce AUROC table, based on all patients (n=188)")}
		library(pROC)
	
		#AUROC SerumIgA	
		roc.result <- roc(IgA.df$BinaryHCC, IgA.df$SerumIgA) #roc(category,predictor)
		auc(roc.result)
		plot(roc.result)
	
		#AUROC SerumIgA50
		roc.result <- roc(IgA.df$BinaryHCC, IgA.df$SerumIgA50) #roc(category,predictor)
		auc(roc.result)
		plot(roc.result)
	
		#manual ROC cutoffs for IgA50
			ROC.df <- IgA.df[order(IgA.df$SerumIgA),][,colnames(IgA.df) %in% c("SerumIgA","SerumIgA50","BinaryHCC")]
			n.events <- sum(ROC.df$BinaryHCC)
	
			for (i in 1:nrow(ROC.df)){
				ROC.df$IgA50HCC[i] <- sum(ROC.df[ROC.df$SerumIgA50 == ROC.df$SerumIgA50[i],]$BinaryHCC)
			}
			
			ROC.df <- ROC.df[!duplicated(as.numeric(ROC.df$SerumIgA50)),][,3:4] #only retain serumIgA50
			
			AUROC.table <- cbind(rbind(rep(0,2),ROC.df),roc.result$sensitivities,roc.result$specificities)
			AUROC.table
	
	
#add serum IgA <=150; >150 & < 400 ; >= 400
	SerumIgA3bin <- rep(0, nrow(IgA.df))
	for(i in seq_along(SerumIgA3bin)){
		if(IgA.df$SerumIgA[i] <= 150){ SerumIgA3bin[i] <-0 }
		if(IgA.df$SerumIgA[i] > 150 & IgA.df$SerumIgA[i] < 400){ SerumIgA3bin[i] <-1 }
		if(IgA.df$SerumIgA[i] > 400){ SerumIgA3bin[i] <-2 }
	}
	
	IgA.df$SerumIgA3bin <- SerumIgA3bin


# Cox regression analysis 
	#### Age coxph
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ Age, data = IgA.df, method = "efron", robust = F)
		summary(coxph.result)
		summary(coxph.result)$logtest
		summary(coxph.result)$logtest[3] %>% as.numeric
		summary(coxph.result)$waldtest
		summary(coxph.result)$waldtest[3] %>% as.numeric
		coxph.result$wald.test #wald test score
		pie.age <- coxph.result$wald.test #wald test score
	
	
	#### SerumIgA50 coxph
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ SerumIgA50, data = IgA.df, method = "efron", robust = F)
		summary(coxph.result)
		summary(coxph.result)$logtest
		summary(coxph.result)$logtest[3] %>% as.numeric
		summary(coxph.result)$waldtest
		summary(coxph.result)$waldtest[3] %>% as.numeric
		coxph.result$wald.test #wald test score

	#### SerumIgA3bin coxph
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ SerumIgA3bin, data = IgA.df, method = "efron", robust = F)
		summary(coxph.result)
		summary(coxph.result)$logtest
		summary(coxph.result)$logtest[3] %>% as.numeric
		summary(coxph.result)$waldtest
		summary(coxph.result)$waldtest[3] %>% as.numeric
		coxph.result$wald.test #wald test score
		pie.serumIgA3bin <- coxph.result$wald.test #wald test score
	
	#### Fibrosis coxph
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ Fibrosis, data = IgA.df, method = "efron", robust = F)
		summary(coxph.result)
		summary(coxph.result)$logtest
		summary(coxph.result)$logtest[3] %>% as.numeric
		summary(coxph.result)$waldtest
		summary(coxph.result)$waldtest[3] %>% as.numeric
		coxph.result$wald.test #wald test score
		pie.fibrosis <- coxph.result$wald.test #wald test score
	
	
	#### HIV
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ HIV, data = IgA.df, method = "efron", robust = F)
		summary(coxph.result)
		summary(coxph.result)$logtest
		summary(coxph.result)$logtest[3] %>% as.numeric
		summary(coxph.result)$waldtest
		summary(coxph.result)$waldtest[3] %>% as.numeric
		coxph.result$wald.test #wald test score
	
	
	#### HCV
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ HCV, data = IgA.df, method = "efron", robust = F)
		summary(coxph.result)
		summary(coxph.result)$logtest
		summary(coxph.result)$logtest[3] %>% as.numeric
		summary(coxph.result)$waldtest
		summary(coxph.result)$waldtest[3] %>% as.numeric
		coxph.result$wald.test #wald test score
	
	#### HBV
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ HBV, data = IgA.df, method = "efron", robust = F)
		summary(coxph.result)
		summary(coxph.result)$logtest
		summary(coxph.result)$logtest[3] %>% as.numeric
		summary(coxph.result)$waldtest
		summary(coxph.result)$waldtest[3] %>% as.numeric
		coxph.result$wald.test #wald test score
	
	#### HBV-HCV
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ HBVHCV, data = IgA.df, method = "efron", robust = F)
		summary(coxph.result)
		summary(coxph.result)$logtest
		summary(coxph.result)$logtest[3] %>% as.numeric
		summary(coxph.result)$waldtest
		summary(coxph.result)$waldtest[3] %>% as.numeric
		coxph.result$wald.test #wald test score
	
	#### HBV-delta
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ HBVDelta, data = IgA.df, method = "efron", robust = F)
		summary(coxph.result)
		summary(coxph.result)$logtest
		summary(coxph.result)$logtest[3] %>% as.numeric
		summary(coxph.result)$waldtest
		summary(coxph.result)$waldtest[3] %>% as.numeric
		coxph.result$wald.test #wald test score
	
	
	#### Gender
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ Gender, data = IgA.df, method = "efron", robust = F)
		summary(coxph.result)
		summary(coxph.result)$logtest
		summary(coxph.result)$logtest[3] %>% as.numeric
		summary(coxph.result)$waldtest
		summary(coxph.result)$waldtest[3] %>% as.numeric
		coxph.result$wald.test #wald test score
	
	
	
	#### SerumIgA3bin +Age + Fibrosis coxph
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ SerumIgA3bin + Age + Fibrosis, data = IgA.df, method = "efron", robust = F)
		summary(coxph.result)
		summary(coxph.result)$logtest
		summary(coxph.result)$logtest[3] %>% as.numeric
		summary(coxph.result)$waldtest
		summary(coxph.result)$waldtest[3] %>% as.numeric
		coxph.result$wald.test #wald test score


	
		#pie graph
			library(ggsci)
			color.vector <- pal_npg("nrc")(3) #requires ggsci package
	
			pie.df <- data.frame(
				group = c("serum IgA (binned <150; 150-400; >=400)", "fibrosis F0/1;F2;F3;F4", "Age"),
				value = c(pie.serumIgA3bin, pie.fibrosis, pie.age)
			)
			
			p=ggplot(pie.df, aes(x="", y=value, fill=group))+
				geom_bar(width = 1, stat = "identity") + 
				coord_polar("y", start=0) +
				theme_void() + scale_fill_manual(values=c(color.vector))
			p

#Akaike Information Criterion
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ Fibrosis, data = IgA.df, method = "efron", robust = F)
		coxph.result ;AIC(coxph.result)
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ Fibrosis + Age  , data = IgA.df, method = "efron", robust = F)
		coxph.result ;AIC(coxph.result)
		coxph.result <- coxph(formula = Surv(CensoredTime, BinaryHCC) ~ Fibrosis + SerumIgA3bin + Age , data = IgA.df, method = "efron", robust = F)
		coxph.result ;AIC(coxph.result)		
		AIC(coxph.result) #indien daalt is het beter bij multivarate cox regression
		
	predict.coxph <- predict(coxph.result)
	IgA.df$predict <- predict.coxph


#Kaplan-Meier curve ASF (Age/Serum IgA/Fibrosis) model: Age >=40 AND Fibrosis = F3/F4 AND serum IgA >150+ mg/dL  is high-risk

	#stratify ASF: Age, Serum IgA and Fibrosis
		ASF <- rep(0,nrow(IgA.df)) 
		for(i in 1:length(ASF)){
			if( IgA.df$Age[i] >= 40 & 	IgA.df$SerumIgA[i] > 150 & IgA.df$F3F4[i] == "Yes") {
			ASF[i] <- "ASFyes"	
			}else{ASF[i] <- "ASFno"}
		}
		IgA.df <- cbind(IgA.df, ASF)

	#create data frame for kaplan meier	  
		km.df <- cbind(IgA.df$BinaryHCC, IgA.df$CensoredTime, IgA.df$ASF)
		colnames(km.df) <- c("HCC","time","ASF")
		km.df <- as.data.frame(km.df)
		
		set.time  = "time"    #column name of time
		set.event = "HCC"   #column name 
	
	#run survfit
		survfit.result      <- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ ASF, data = km.df)
	
	#create truncated Kaplan Meier curve when patients at risk drops below 10% for visualization purpose (not for p-value!)
		#store p-values
			pval.event.plot <- surv_pvalue(survfit.result)$pval

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
		
		#visualize			
			p=ggsurvplot(survfit.result,
				fun= "event",
				pval = paste0("high risk VS low risk \n p=", signif(surv_pvalue(survfit.result)$pval,3)),
				pval.coord = c(0,0.75),
				pval.size=7,
				conf.int = FALSE,
				risk.table = TRUE, risk.table.fontsize=6,
				risk.table.height=0.3,
				risk.table.y.text = TRUE,
				font.x=20,
				font.y=20,
				font.tickslab=20,
				palette = c("green2","orange"),
				xlim = c(0,cutoff.time.x.axis), 
				break.time.by=3,
				legend = 'none',
				xlab="Time (years)",
				ylab="HCC incidence"
				)
				p$table <- p$table + 
						theme(axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20)) +
						theme(axis.title.y = element_text(size = 20)) + ylab("HCC risk") + 
						scale_y_discrete(labels=c('high', 'low'))
			p


		

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

	#create data frame for Kaplan-Meier	  
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
		
		#visualise kaplan meier plot
			p=ggsurvplot(survfit.result,
				fun= "event",
				pval = paste0("high VS low p=",pval.low.vs.high,"\n", "high VS int. p=",pval.int.vs.high,"\n" , "int. VS low p=",pval.low.vs.int ),   
				pval.coord = c(0,0.90), 
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
			p

	
#KM curve AS2F model: Age >=40 AND F3/F4 AND >=400 mg/dL serum IgA is high-risk, >=40 Age AND F3/F4 AND 150<IgA<400 is intermediate-risk
	#add stratification to data
		IgA.df$ASF <- NULL
		ASF <- rep(0,nrow(IgA.df)) 
		for(i in 1:length(ASF)){
			if(IgA.df$Age[i] >= 40 & IgA.df$SerumIgA[i] >= 400 & IgA.df$F3F4[i] == "Yes") { ASF[i] <- "ASFpoor"	
			}else{
			if(IgA.df$Age[i] >= 40 & IgA.df$SerumIgA[i] > 150 & IgA.df$SerumIgA[i] < 400 & IgA.df$F3F4[i] == "Yes") { ASF[i] <- "ASFintermediate"	
			}else{
			ASF[i] <- "ASFgood"}
			}
		}
		IgA.df <- cbind(IgA.df, ASF)

	#create data frame for Kaplan-Meier	 	  
		km.df <- cbind(IgA.df$BinaryHCC, IgA.df$CensoredTime, IgA.df$ASF)
		colnames(km.df) <- c("HCC","time","ASF")
		km.df <- as.data.frame(km.df)
	
		set.time  = "time"    #column name of time
		set.event = "HCC"   #column name 
		
	#run survfit
		survfit.result       				<- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ ASF, data = km.df)
		survfit.intermediate.removed.result	<- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ ASF, data = km.df[which(km.df$ASF != "ASFintermediate"),])
		survfit.low.removed.result			<- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ ASF, data = km.df[which(km.df$ASF != "ASFgood"),])
		survfit.high.removed.result			<- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ ASF, data = km.df[which(km.df$ASF != "ASFpoor"),])

	#create truncated Kaplan Meier curve when patients at risk drops below 10% for visualization purpose (not for p-value!)
		#store p-values
			p.pairwise 			<- pairwise_survdiff(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ ASF, data = km.df, p.adjust.method = "none")
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
		
		#visualize Kaplan Meier
				p=ggsurvplot(survfit.result,
					fun= "event",
					pval = paste0("high VS low p=",pval.low.vs.high,"\n", "high VS int. p=",pval.int.vs.high,"\n" , "int. VS low p=",pval.low.vs.int ),   
					pval.coord = c(0,0.70), 
					pval.size=7,
					conf.int = FALSE,
					risk.table = TRUE, risk.table.fontsize=6,
					risk.table.height=0.3,
					risk.table.y.text = TRUE,
					font.x=20,
					font.y=20,
					font.tickslab=20,
					palette = c("green2", "darkgrey","purple"),
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
							scale_y_discrete(labels=c('high', 'int.', 'low'))
				p


