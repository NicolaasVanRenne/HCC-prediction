################################################################
#
# Calculate 557 prognostic gene signature 
#
#################################################################
#note: The gene signature can be calculated using the SurvivalGene and LoocvSurvival algorithms, created by Yujin Hoshida.
#note: The algorithms have been described earlier (Hoshida et al. 2008, NEJM) and were until recently freely avalaible on genepattern archive (http://gparc.org/).
#note: Documentation of the algorithms is still available on GenePattern GitHub https://github.com/genepattern/gparc-module-docs
#note: Both algorithms were integrated here into a single pipeline.

#note: If you use these algorithms, please cite "Hoshida et al. 2008, NEJM, PMID:18923165"

#note: Since these algorithms are redistributed from GenePattern, the original copyright notice is included here:

#	----------------------------------------------------------------
#	       *** GENEPATTERN LICENSE AGREEMENT ***
#	                   
#	----------------------------------------------------------------
#	GenePattern is distributed under the following BSD-style license:
#	                            
#	Copyright (c) 2003-2022 Regents of the University of California and Broad Institute. All rights reserved.
#	                            
#	Redistribution and use in source and binary forms, with or without modification, 
#	are permitted provided that the following conditions are met:
#	                            
#	     1. Redistributions of source code must retain the above copyright notice, 
#	        this list of conditions and the following disclaimer.
#	                            
#	     2. Redistributions in binary form must reproduce the above copyright notice, 
#	        this list of conditions and the following disclaimer in the documentation 
#	        and/or other materials provided with the distribution.
#	                            
#	     3. Neither the names of the Broad Institute, Inc., Massachusetts Institute of Technology 
#	        nor the names of its contributors may be used to endorse or promote products derived 
#	        from this software without specific prior written permission.
#	                            
#	THIS SOFTWARE IS PROVIDED AS IS.  BROAD AND MIT MAKE NO EXPRESS OR IMPLIED 
#	REPRESENTATIONS OR WARRANTIES OF ANY KIND REGARDING THE SOFTWARE AND COPYRIGHT, 
#	INCLUDING, BUT NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#	PARTICULAR PURPOSE, CONFORMITY WITH ANY DOCUMENTATION, NONINFRINGEMENT, OR THE 
#	ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. IN NO EVENT SHALL
#	BROAD, MIT, THE COPYRIGHT HOLDERS, OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
#	INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
#	BUT NOT LIMITED TO PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
#	OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
#	WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF, HAVE REASON 
#	TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF SUCH DAMAGE.
#	                            
#	If, by operation of law or otherwise, any of the aforementioned warranty disclaimers 
#	are determined inapplicable, your sole remedy, regardless of the form of action, 
#	including, but not limited to, negligence and strict liability, shall be replacement
#	of the software with an updated version if one exists.


##############################################################################################
# step 1. SurvivalGene
##############################################################################################

# note: the gene expression data should be loaded as a .gct file (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats)
# add this on line 1 of the gct file: '#1.2'
# add this on line 2 of the gct file: '20536	56'


#SETTINGS
	setwd("C:/my_dir")
	
	input.filename.gct="data_files/collapsed_data/GSE237330_training_RPM.gct" 	#gene expression file in gct format
	input.filename.clinical="data_files/meta_data/training_clinical_HCC.txt" # .txt file with time and events
	output.file = "training_RPM_filtered_signature"									#output file name
	
	set.p.value = 0.1 #p.value cutoff for LOOCV
	
	#set random seed version to old (R version < 3.6) or new (R version > 3.6)
		choose.seed = "old" #"old" or "new": set random seed version to old (R version < 3.6; reproduces GenePattern) or new (R version > 3.6)
		if(choose.seed == "old"){RNGkind(sample.kind = "Rounding")}else{} #"old" necessary to reproduce GenePattern results
		
	
#RUN ALGORITHM
	
	##############################################################
	# SurvivalGene.R (ver.2)                  Aug.13,2008
	#    (modified for GenePattern Oct 13, 2011)
	#     
	#     input:  (1) gct
	#             (2) clinical data
	#                *.txt, 1st column is sample name in the same
	#                 order with .gct
	#                *
	#     output: (1) result
	#             (2) mean & sd of empirical dist for the statistic
	#                 (used for loocv.nn.R & nn.surv.iter.Kfold.R)
	#                 Shapiro-Wilk normality test (W & p-value)
	##############################################################
	
	{
	time.field="time"
	censor.field="status"
	
	statistic.selection="cox.score" #  "cox.score", "coxph"
	trim.percent.2.side=0.00
	nperm=1000
	rnd.seed=56438219
	
	emp.stat.dist="T" 
	
	

	
	{
	
	# excluded from arguments for GenePattern
	check.names.gct="F"
	
	# for GenePattern
	
	trim.percent.2.side<-as.numeric(trim.percent.2.side)
	rnd.seed <- as.numeric(rnd.seed)
	nperm <- as.numeric(nperm)
	
	set.seed(rnd.seed)
	
	# read input data
	
	if (check.names.gct=="F"){
		exp.all <- read.delim(input.filename.gct,header=T,skip=2,check.names=F)
	}else{
		exp.all <- read.delim(input.filename.gct,header=T,skip=2)
	}
	
	## file format check
	if (regexpr(".gct$",input.filename.gct)==-1){
		stop("### Gene expression data should be .gct format! ###")
	}
	
	probeid <- as.vector(exp.all[,1])
	gene.names <- as.vector(exp.all[,2])
	num.samples<-(length(exp.all[1,])-2)
	num.genes <- length(exp.all[,1])
	exp.all<-exp.all[-c(1:2)]
	sample.names.all <- colnames(exp.all)
	rownames(exp.all) <- probeid
	
	clin.all <- read.delim(input.filename.clinical,header=T)
	
	eval(parse(text=paste("time <- clin.all$",time.field,sep="")))
	eval(parse(text=paste("censor <- clin.all$",censor.field,sep="")))
	
	sample.names.clin <- as.vector(clin.all[,1])
	
	## check sample order btwn exp & clin
	for (i in 1:num.samples){
		if (sample.names.all[i]!=sample.names.clin[i]){
		stop("### Sample names don't match! ###")
		}
	}
	
	# trim survival time outliers
	
	num.trim.1.side <- round((num.samples*trim.percent.2.side+0.5)/2)
	order.survival <- order(time)
	trimmed.ordered.index <- order.survival[(num.trim.1.side+1):(num.samples-num.trim.1.side)]
	trimmed.time <- time[trimmed.ordered.index]
	trimmed.censor <- censor[trimmed.ordered.index]
	trimmed.exp <- exp.all[,trimmed.ordered.index]
	num.samples.trimmed <- length(trimmed.time)
	
	# set permutation matrices
	
	perm.stat <- matrix(0,nrow=num.genes,ncol=(nperm+1))
	perm.index <- matrix(0,nrow=num.samples.trimmed,ncol=nperm)
	perm.index.vector <- seq(1:num.samples.trimmed)
	for (i in 1:nperm){
		perm.index[,i] <- sample(perm.index.vector)
	}
	
	# cases at risk & death
	event.index <- which(trimmed.censor==1)
	event.time <- unique(trimmed.time[event.index])
	num.event.time <- length(event.time)
	
	# compute statistic
	if (statistic.selection=="cox.score"){
	
		# original statistics
		numerator <- 0
		denominator <- 0
	
		for (t in 1:num.event.time){
		index.at.risk <- which(trimmed.time>=event.time[t])
		num.at.risk <- length(index.at.risk)
		index.death <- which(trimmed.time==event.time[t])
		num.death <- length(index.death)
	
		if (is.null(dim(trimmed.exp[,index.at.risk]))==T){
			if (num.death==1){
			numerator <- numerator + trimmed.exp[,index.death]-num.death*trimmed.exp[,index.at.risk]
			}
		# if (num.death>1){
		#   numerator <- numerator + apply(trimmed.exp[,index.death],1,sum)-num.death*apply(trimmed.exp[,index.at.risk],1,mean)
		# }
			denominator <- denominator + (num.death/num.at.risk)*as.vector((trimmed.exp[,index.at.risk]-as.vector(trimmed.exp[,index.at.risk]))^2)
		}else{
			if (num.death==1){
			numerator <- numerator + trimmed.exp[,index.death]-num.death*apply(trimmed.exp[,index.at.risk],1,mean)
			}
			if (num.death>1){
			numerator <- numerator + apply(trimmed.exp[,index.death],1,sum)-num.death*apply(trimmed.exp[,index.at.risk],1,mean)
			}
			denominator <- denominator + (num.death/num.at.risk)*apply((trimmed.exp[,index.at.risk]-apply(trimmed.exp[,index.at.risk],1,mean))^2,1,sum)
	#        denominator <- denominator + (num.death/num.at.risk)*as.vector(apply((trimmed.exp[,index.at.risk]-as.vector(apply(trimmed.exp[,index.at.risk],1,mean)))^2,1,sum))
		}      
		}
		perm.stat[,1] <-numerator/sqrt(denominator)
	
		# permutation startistic
	
		for (p in 1:nperm){
	
		print(paste("# Start permutation ",p,sep=""))
	
		perm.trimmed.exp <- trimmed.exp[,perm.index[,p]]
	
		numerator <- 0
		denominator <- 0
	
		for (t in 1:num.event.time){
			index.at.risk <- which(trimmed.time>=event.time[t])
			num.at.risk <- length(index.at.risk)
			index.death <- which(trimmed.time==event.time[t])
			num.death <- length(index.death)
	
			if (is.null(dim(perm.trimmed.exp[,index.at.risk]))==T){
			if (num.death==1){
				numerator <- numerator + perm.trimmed.exp[,index.death]-num.death*perm.trimmed.exp[,index.at.risk]
			}
			# if (num.death>1){
			#   numerator <- numerator + apply(perm.trimmed.exp[,index.death],1,sum)-num.death*apply(perm.trimmed.exp[,index.at.risk],1,mean)
			# }
			denominator <- denominator + (num.death/num.at.risk)*(perm.trimmed.exp[,index.at.risk]-perm.trimmed.exp[,index.at.risk])^2
			}else{
			if (num.death==1){
				numerator <- numerator + perm.trimmed.exp[,index.death]-num.death*apply(perm.trimmed.exp[,index.at.risk],1,mean)
			}
			if (num.death>1){
				numerator <- numerator + apply(perm.trimmed.exp[,index.death],1,sum)-num.death*apply(perm.trimmed.exp[,index.at.risk],1,mean)
			}
			denominator <- denominator + (num.death/num.at.risk)*apply((perm.trimmed.exp[,index.at.risk]-apply(perm.trimmed.exp[,index.at.risk],1,mean))^2,1,sum)
			}
		}
		perm.stat[,(p+1)] <-numerator/sqrt(denominator)
		}   # permutation statistics END
	
		# p-value
	
		rank.perm.stat <- nperm+1-t(apply(abs(perm.stat),1,rank,ties.method="min"))
		nominal.p <- rank.perm.stat[,1]/nperm
	
		rank.nom.p <- rank(nominal.p,ties.method="min")
		BH.FDR <- nominal.p*num.genes/rank.nom.p
		BH.FDR[BH.FDR>1] <- 1
		FWER <- nominal.p*num.genes
		FWER[FWER>1] <- 1
	
	} # cox.score END
	
	if (statistic.selection=="coxph"){
		library(survival)
		statistic <- p.value <- vector(length=num.genes,mode="numeric")
		for (g in 1:num.genes){
		print(paste("# Start gene ",g,"/",num.genes,sep=""))
		coxfit <- coxph(Surv(trimmed.time,trimmed.censor)~as.vector(as.matrix(trimmed.exp[g,])))
		statistic[g] <- as.numeric(coxfit$coef)
		p.value[g] <- as.numeric(unlist(summary(coxfit))$coefficients5)
		}
	
		rank.nom.p <- rank(p.value,ties.method="min")
		BH.FDR <- p.value*num.genes/rank.nom.p
		BH.FDR[BH.FDR>1] <- 1
		FWER <- p.value*num.genes
		FWER[FWER>1] <- 1
	
	}  # coxph END
	
	# output
	if (statistic.selection=="cox.score"){
		statistic <- perm.stat[,1]
		output <- cbind(probeid,gene.names,statistic,nominal.p,BH.FDR,FWER)
	}
	if (statistic.selection=="coxph"){
		output <- cbind(probeid,gene.names,statistic,p.value,BH.FDR,FWER)
	}
	write.table(output,paste(output.file,".txt",sep=""),quote=F,sep="\t",row.names=F)
	
	if (emp.stat.dist=="T" & statistic.selection=="cox.score"){
		mean <- apply(perm.stat[,2:(nperm+1)],1,mean)
		sd <- apply(perm.stat[,2:(nperm+1)],1,sd)
	#    mean <- apply(perm.stat[,2:(nperm+1)],1,mean,na.rm=T)
	#    sd <- apply(perm.stat[,2:(nperm+1)],1,sd,na.rm=T)
		sk.test.W <- sk.test.p <- vector(length=num.genes,mode="numeric")
		for (i in 1:num.genes){
		sk.test <- shapiro.test(as.vector(perm.stat[i,2:(nperm+1)]))
		sk.test.W[i] <- sk.test$statistic
		sk.test.p[i] <- sk.test$p.value
		}
	
		emp.stat.output <- cbind(probeid,mean,sd,sk.test.W,sk.test.p)
		write.table(emp.stat.output,paste(output.file,"_emp.stat.txt",sep=""),quote=F,sep="\t",row.names=F)
	}
	
	#  return(list(output,perm.stat))
	
	} # main END
}

	#############################################################
	# loocv.nn_cox_GP.R (ver.2.1)                Feb.12, 2009
	#    (modified for GenePattern   Oct 15, 2011)
	#
	#   LOOCV using Nearest Template Prediction (NTP)
	#
	#   Gene ranking metric:
	#     Survival cls (only NTP): Cox score, Cox regression, samr()
	#
	#   Input: gct (1st column: probeid, 2nd column: symbol)
	#               cls (original cls defined w/o validation, should start from "1")
	#               clinical data (txt w/ header)
	#   Output: Usage of up & down genes
	#              LOOCV prediction result
	#                  (.txt, .cls. dChip sample info)
	#############################################################

{	
	input.dataset=input.filename.gct
	input.clinical.data=input.filename.clinical
	
	# for CoxScore
	time.field="time"
	censor.field="status"
	trim.percent.2.side=0
	cox.score.sig= set.p.value  # 2-sided p-value. >=1 is regarded as # of marker genes
						# assume normal dist
	emp.cox.file=paste(output.file,"_emp.stat.txt",sep="")   # emripical dist of Cox score to standardize
	
	# Prediction
	temp.nn.wt="T" # "T" reproduces Hoshida 2008 NEJM
	dist.selection="cosine"  # "correlation" or "cosine"
	
	# number of resampling to generate null dist for prediction stat
	nresample=1000
	
	# random seed
	rnd.seed=4675921
	
	output.file.2=paste0(output.file,"_LOOCV_Survival_NTP")
	
	
	{
		# Auxilary functions
		
		read.gct<-function(filename="NULL")
		{
		
		if (regexpr(".gct$",filename)==-1){
			stop("### input data should be .gct file! ###")
		}
		
		data<-read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T,check.names=F)
		data<-data[-1]
		
		return(data)
		}
		
		read.gct.gene.names<-function(filename="NULL")
		{
		data<-read.delim(filename, header=T, skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T)
		gene<-as.vector(data[,1])
		
		return(gene)
		}
		
		read.cls<-function(filename="NULL")
		{
		cls<-as.vector(as.matrix(read.delim(filename,header=F,sep=" ",skip=2)))
		return(cls)
		}
		
		cosine<-function(v1,v2)
		{
		cosine.dist<-sum(v1*v2)/sqrt(sum(v1^2)*sum(v2^2))
		return(cosine.dist)
		}
		
		ttest <- function(cls.1.exp,cls.2.exp,ttest.type="welch",sigma.correction="GeneCluster",
						num.cls.1=num.cls.1,num.cls.2=num.cls.2)
		{
		cls.1.mean <- apply(cls.1.exp,1,mean)
		cls.2.mean <- apply(cls.2.exp,1,mean)
		cls.1.sd <- apply(cls.1.exp,1,sd)
		cls.2.sd <- apply(cls.2.exp,1,sd)
		
		# GeneCluster type small sigma fix
		
		if (sigma.correction == "GeneCluster") {
			cls.1.sd <- ifelse(0.2*abs(cls.1.mean) < cls.1.sd, cls.1.sd, 0.2*abs(cls.1.mean))
			cls.1.sd <- ifelse(cls.1.sd == 0, 0.2, cls.1.sd)
			cls.2.sd <- ifelse(0.2*abs(cls.2.mean) < cls.2.sd, cls.2.sd, 0.2*abs(cls.2.mean))
			cls.2.sd <- ifelse(cls.2.sd == 0, 0.2, cls.2.sd)
			gc()
		}
		
		ttest <- (cls.1.mean - cls.2.mean)/sqrt(((cls.1.sd^2)/num.cls.1) + ((cls.2.sd^2)/num.cls.2))
		return(ttest)
		}
		
		CoxScore <- function(loo.dataset,loo.time,loo.censor,trim.percent.2.side=0)
		{
		
		#  probeid <- as.vector(exp.all[,1])
		#  gene.names <- as.vector(exp.all[,2])
		num.samples<-length(loo.dataset[1,])
		num.genes <- length(loo.dataset[,1])
		#  sample.names.all <- colnames(exp.all)
		#  rownames(exp.all) <- probeid
		
		## check sample order btwn exp & clin
		#  for (i in 1:num.samples){
		#    if (sample.names.all[i]!=sample.names.clin[i]){
		#      stop("### Sample names don't match! ###")
		#    }
		#  }
		
		# trim survival time outliers
		
		num.trim.1.side <- round((num.samples*trim.percent.2.side+0.5)/2)
		order.survival <- order(loo.time)
		
		trimmed.ordered.index <- order.survival[(num.trim.1.side+1):(num.samples-num.trim.1.side)]
		
		trimmed.time <- loo.time[trimmed.ordered.index]
		trimmed.censor <- loo.censor[trimmed.ordered.index]
		trimmed.exp <- loo.dataset[,trimmed.ordered.index]
		num.samples.trimmed <- length(trimmed.time)
		
		# cases at risk & death
		event.index <- which(trimmed.censor==1)
		event.time <- unique(trimmed.time[event.index])
		num.event.time <- length(event.time)
		
		# compute statistic
		numerator <- 0
		denominator <- 0
		
		for (t in 1:num.event.time){
			index.at.risk <- which(trimmed.time>=event.time[t])
			num.at.risk <- length(index.at.risk)
			index.death <- which(trimmed.time==event.time[t])
			num.death <- length(index.death)
		
			if (is.null(dim(trimmed.exp[,index.at.risk]))){
			if (num.death==1){
			numerator <- numerator + trimmed.exp[,index.death]-num.death*trimmed.exp[,index.at.risk]
			}
			# if (num.death>1){
			#   numerator <- numerator + apply(trimmed.exp[,index.death],1,sum)-num.death*apply(trimmed.exp[,index.at.risk],1,mean)
			# }
			denominator <- denominator + (num.death/num.at.risk)*as.vector((trimmed.exp[,index.at.risk]-as.vector(trimmed.exp[,index.at.risk]))^2)
			}else{
			if (num.death==1){
			numerator <- numerator + trimmed.exp[,index.death]-num.death*apply(trimmed.exp[,index.at.risk],1,mean)
			}
			if (num.death>1){
				numerator <- numerator + apply(trimmed.exp[,index.death],1,sum)-num.death*apply(trimmed.exp[,index.at.risk],1,mean)
			}
			denominator <- denominator + (num.death/num.at.risk)*apply((trimmed.exp[,index.at.risk]-apply(trimmed.exp[,index.at.risk],1,mean))^2,1,sum)
			}
		}
		output <-numerator/sqrt(denominator)
		return(output)
		}
		
		# MAIN PART
		
		
		# file I/O
		input.cls=NULL
		cls.label=NULL
		
		# marker selection: removed for GP version
		marker.selection="cox.score"  # "ttest","snr","cox.score","samr.surv" (snr not implemented yet): allow just t-test & Cox score for GP version
		
		## t-test-related arguments: removed for GP version
		t.stat.sig=0.05  # 2-sided
		sigma.correction="none" # "GeneCluster"
		
		## for samr(): removed from arguments for GP version
		nperm.samr=100
		delta=.4
		
		# prediction-related arguments: removed from arguments for GP version
		pred.method="temp.nn"     # "temp.nn" or "knn"
		within.sig="F"
		histgram.cosine="F"
		within.sig="F"
		histgram.cosine="F"
		
		# for GenePattern
		t.stat.sig <- as.numeric(t.stat.sig)
		trim.percent.2.side <- as.numeric(trim.percent.2.side)
		cox.score.sig <- as.numeric(cox.score.sig)
		nperm.samr <- as.numeric(nperm.samr)
		delta <- as.numeric(delta)
		nresample <- as.numeric(nresample)
		rnd.seed <- as.numeric(rnd.seed)
		
		set.seed(rnd.seed)          # set random seed
		
		# packages
		
		if (marker.selection=="samr.surv"){
			require(samr,quietly=T)
		}
		
		# read expression data
		
		dataset<-read.gct(input.dataset)
		dataset.row.norm<-(dataset-apply(dataset,1,mean,na.rm=T))/apply(dataset,1,sd,na.rm=T)
		#print(dataset.row.norm[1:3,1:3])
		gene.names<-read.gct.gene.names(input.dataset)
		probeid<-rownames(dataset)
		num.gene<-length(dataset[,1])
		gene.id<-rownames(dataset)
		num.sample<-length(dataset[1,])
		sample.names<-colnames(dataset)
		
		# read cls clinical data
		
		if (is.null(cls.label) & is.null(input.cls) & (marker.selection=="ttest" | marker.selection=="snr")){
			stop("# Please provide cls information!!!")
		}
		if (is.null(cls.label) & (marker.selection=="ttest" | marker.selection=="snr")){
			cls<-read.cls(input.cls)
			if(length(cls)!=num.sample){
			stop("# of samples in .gct and .cls don't match!!!")
			}
			cls.label.1 <- as.vector(unique(cls)[1])
			cls.label.2 <- as.vector(unique(cls)[2])
		}
		if (is.null(input.cls) & (marker.selection=="ttest" | marker.selection=="snr")){
			clin<-read.delim(file=input.clinical.data,header=T)
			eval(parse(text=paste("cls<-clin$",cls.label,sep="")))
			if(length(cls)!=num.sample){
			stop("# of samples in .gct and cls don't match!!!")
			}
			eval(parse(text=paste("cls.labels <- unique(clin$",cls.label,")",sep="")))
			cls.label.1 <- as.vector(cls.labels[1])
			cls.label.2 <- as.vector(cls.labels[2])
		}
		
		if (marker.selection=="cox.score" | marker.selection=="samr.surv"){
			clin<-read.delim(file=input.clinical.data,header=T)
			time<-eval(parse(text=paste("clin$",time.field,sep="")))
			censor<-eval(parse(text=paste("clin$",censor.field,sep="")))
			if (is.null(time) | is.null(censor)){
			stop("Please provide survival information!")
			}
		}
		
		# read empirical Cox score (if specified)
		if (emp.cox.file=="NULL"){
			emp.cox.file <- NULL
		}
		if (is.null(emp.cox.file)){
		}else{
			emp.cox <- read.delim(file=emp.cox.file,header=T)
			for (g in 1:num.gene){
			if (probeid[g]!=emp.cox[g,1]){
				stop("### ProbeID of empirical Cox score doesn't match! ###")
			}
			}
		}
		
		# output objects
		
		if (marker.selection=="ttest" | marker.selection=="snr"){
			pred.result<-as.data.frame(matrix(0,nrow=num.sample,ncol=7))
			pred.result[,1]<-sample.names
			pred.result[,2]<-cls
			colnames(pred.result)<-c("SampleName","OriginalCls","PredictedCls","CosineDist","nominal.p","BH.FDR","Bonferroni")
		
			up.genes.feature.table<-cbind(probeid,gene.names)
			colnames(up.genes.feature.table)<-c("probeid","GeneName")
		
			down.genes.feature.table<-cbind(probeid,gene.names)
			colnames(down.genes.feature.table)<-c("probeid","GeneName")
		}
		if (marker.selection=="cox.score" | marker.selection=="samr.surv"){
			pred.result<-as.data.frame(matrix(0,nrow=num.sample,ncol=7))
			pred.result[,1]<-sample.names
			colnames(pred.result)<-c("SampleName","OriginalCls","PredictedCls","CosineDist","nominal.p","BH.FDR","FWER")
		
			poor.genes.feature.table<-cbind(probeid,gene.names)
			colnames(poor.genes.feature.table)<-c("probeid","GeneName")
		
			good.genes.feature.table<-cbind(probeid,gene.names)
			colnames(good.genes.feature.table)<-c("probeid","GeneName")
		#print("pred.result is set")
		}
		
		resample.stat <- matrix(0,nrow=num.sample,ncol=nresample*2)
		
		# threshold for marker gene selection
		
		if (t.stat.sig<1){
			t.stat.threshold <- qt((1-(t.stat.sig/2)),(num.sample-1))   # t-stat significance threshold
		}else{
			num.markers <- t.stat.sig
		}
		if (cox.score.sig<1){
			cox.score.threshold <- qnorm((1-(cox.score.sig/2)))   # cox score significance threshold
		}else{
			num.markers <- cox.score.sig
		}
		
		# LOOCV loop
		
		for (loo in 1:num.sample){
		
			print(paste("LOOCV sample ",loo,"/",num.sample,sep=""))
		
			# generate LOO dataset/cls and sample/cls
			loo.dataset<-dataset[-loo]
			loo.sample.data<-dataset.row.norm[,loo]
			loo.sample.data<-cbind(probeid,loo.sample.data)
		
			if (marker.selection=="ttest" | marker.selection=="snr"){
			loo.cls<-cls[-loo]
			loo.sample.cls<-cls[loo]
			}
			if (marker.selection=="cox.score" | marker.selection=="samr.surv"){
			loo.time<-time[-loo]
			loo.censor<-censor[-loo]
			}
		
			# gene selection
		
			if (marker.selection=="cox.score"){    # marker selection: Cox score START
		
			# rank genes & pick significant genes
			gene.rank.metric.vector <- CoxScore(loo.dataset,loo.time,loo.censor,trim.percent.2.side=trim.percent.2.side)
		#print(paste("grmv: ",gene.rank.metric.vector[1],sep=""))
			# normalize Cox score
			if (is.null(emp.cox.file)){
			}else{
				gene.rank.metric.vector <- (gene.rank.metric.vector-as.vector(emp.cox[,2]))/as.vector(emp.cox[,3])
			}
		#print(paste("mean: ",emp.cox[1,2],sep=""))
		#print(paste("sd: ",emp.cox[1,3],sep=""))
		#print(gene.rank.metric.vector[1])
		
			# select markers
			if (cox.score.sig < 1){    # in case of significance threshold
				index.poor.markers <- which(gene.rank.metric.vector > cox.score.threshold)
				index.good.markers <- which(gene.rank.metric.vector < -cox.score.threshold)
			}
		
			if (cox.score.sig >= 1){     # fixed number marker genes
				order.probeid <- order(gene.rank.metric.vector,decreasing=T)
				index.poor.markers <- order.probeid[1:num.markers]
				index.good.markers <- order.probeid[(num.gene-num.markers+1):num.gene]
			}
		
			probeid.poor.markers <- probeid[index.poor.markers]
			probeid.good.markers <- probeid[index.good.markers]
			marker.stat <- gene.rank.metric.vector[c(index.poor.markers,index.good.markers)]
		
		print(paste(length(index.poor.markers)," poor marker genes",sep=""))
		print(paste(length(index.good.markers)," good marker genes",sep=""))
		
			# make each part of template (i.e. "poor" & "good")
			if (length(index.poor.markers)==0){
		
				eval(parse(text=paste(sample.names[loo],"<-NA",sep="")))
				poor.genes.feature.table<-eval(parse(text=paste("cbind(poor.genes.feature.table,",sample.names[loo],")",sep="")))
		
				poor.genes.temp<-NULL
			}else{
				poor.genes.id<-as.data.frame(cbind(probeid.poor.markers,"poor"))
				colnames(poor.genes.id)<-c("probeid",sample.names[loo])
				poor.genes.feature.table<-merge(poor.genes.feature.table,poor.genes.id,all.x=T)
		
				poor.genes.temp<-cbind(probeid.poor.markers,1,0)
			}
		
			if (length(index.good.markers)==0){
		
				eval(parse(text=paste(sample.names[loo],"<-NA",sep="")))
				good.genes.feature.table<-eval(parse(text=paste("cbind(good.genes.feature.table,",sample.names[loo],")",sep="")))
		
				good.genes.temp<-NULL
			}else{
				good.genes.id<-as.data.frame(cbind(probeid.good.markers,"good"))
				colnames(good.genes.id)<-c("probeid",sample.names[loo])
				good.genes.feature.table<-merge(good.genes.feature.table,good.genes.id,all.x=T)
		
				good.genes.temp<-cbind(probeid.good.markers,0,1)
			}
		
			} # marker selection: Cox Score END
		
			if (marker.selection=="samr.surv"){           # marker selection: samr START
			data.samr.feed=list(x=loo.dataset,y=loo.time,censoring.status=loo.censor,geneid=gene.names,genenames=probeid)
			samr.obj<-samr(data.samr.feed,resp.type="Survival", nperms=nperm.samr,random.seed=rnd.seed)
		
			delta.table <- samr.compute.delta.table(samr.obj)
		
			siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, data.samr.feed, delta.table)
		
			if (siggenes.table$ngenes.up==0){
		
				eval(parse(text=paste(sample.names[loo],"<-NA",sep="")))
				poor.genes.feature.table<-eval(parse(text=paste("cbind(poor.genes.feature.table,",sample.names[loo],")",sep="")))
		
				poor.genes.temp<-NULL
			}else{
				poor.genes.table<-siggenes.table$genes.up
				poor.genes.id<-as.data.frame(cbind(poor.genes.table[,2],"poor"))
				colnames(poor.genes.id)<-c("probeid",sample.names[loo])
				poor.genes.feature.table<-merge(poor.genes.feature.table,poor.genes.id,all.x=T)
		
				poor.genes.temp<-cbind(poor.genes.table[,2],1,0)
			}
		
			if (siggenes.table$ngenes.lo==0){
		
				eval(parse(text=paste(sample.names[loo],"<-NA",sep="")))
				good.genes.feature.table<-eval(parse(text=paste("cbind(good.genes.feature.table,",sample.names[loo],")",sep="")))
		
				good.genes.temp<-NULL
			}else{
				good.genes.table<-siggenes.table$genes.lo
				good.genes.id<-as.data.frame(cbind(good.genes.table[,2],"good"))
				colnames(good.genes.id)<-c("probeid",sample.names[loo])
				good.genes.feature.table<-merge(good.genes.feature.table,good.genes.id,all.x=T)
		
				good.genes.temp<-cbind(good.genes.table[,2],0,1)
			}
			}   # marker selection: samr.surv END
		
			if (marker.selection=="ttest"){     # marker selection: t-test START
			# make 0/1 cls & exp data for 2-cls comparison
			cls.bin.1 <- cls.bin.2 <- as.vector(loo.cls)
		#print(cls.bin.1)
		#print(cls.label.1)
			cls.bin.1[cls.bin.1==cls.label.1] <- 1
			cls.bin.1[cls.bin.1==cls.label.2] <- 0
			cls.bin.1 <- as.numeric(cls.bin.1)
			index.cls.bin.1 <- which(cls.bin.1==1)
		
			cls.bin.2[cls.bin.2==cls.label.1] <- 0
			cls.bin.2[cls.bin.2==cls.label.2] <- 1
			cls.bin.2 <- as.numeric(cls.bin.2)
			index.cls.bin.2 <- which(cls.bin.2==1)
		
			cls.bin.1.exp <- loo.dataset[,index.cls.bin.1]
			cls.bin.2.exp <- loo.dataset[,index.cls.bin.2]
		
			num.cls.1 <- sum(cls.bin.1)
			num.cls.2 <- sum(cls.bin.2)
		
			# rank genes & pick significant genes
			gene.rank.metric.vector <- ttest(cls.bin.1.exp,cls.bin.2.exp,ttest.type="welch",
												sigma.correction=sigma.correction,
												num.cls.1=num.cls.1,num.cls.2=num.cls.2)
		
			index.up.markers <- which(gene.rank.metric.vector > t.stat.threshold)
			index.down.markers <- which(gene.rank.metric.vector < -t.stat.threshold)
		
			probeid.up.markers <- probeid[index.up.markers]
			probeid.down.markers <- probeid[index.down.markers]
		
			marker.stat <- gene.rank.metric.vector[c(index.up.markers,index.down.markers)]
		
		print(paste(length(index.up.markers)," up genes",sep=""))
		print(paste(length(index.down.markers)," down genes",sep=""))
		
			# make each part of template (i.e. "up" & "down")
			if (length(index.up.markers)==0){
		
				eval(parse(text=paste(sample.names[loo],"<-NA",sep="")))
				up.genes.feature.table<-eval(parse(text=paste("cbind(up.genes.feature.table,",sample.names[loo],")",sep="")))
		
				up.genes.temp<-NULL
			}else{
				up.genes.id<-as.data.frame(cbind(probeid.up.markers,"up"))
				colnames(up.genes.id)<-c("probeid",sample.names[loo])
				up.genes.feature.table<-merge(up.genes.feature.table,up.genes.id,all.x=T)
		
				up.genes.temp<-cbind(probeid.up.markers,1,0)
			}
		
			if (length(index.down.markers)==0){
		
				eval(parse(text=paste(sample.names[loo],"<-NA",sep="")))
				down.genes.feature.table<-eval(parse(text=paste("cbind(down.genes.feature.table,",sample.names[loo],")",sep="")))
		
				down.genes.temp<-NULL
			}else{
				down.genes.id<-as.data.frame(cbind(probeid.down.markers,"down"))
				colnames(down.genes.id)<-c("probeid",sample.names[loo])
				down.genes.feature.table<-merge(down.genes.feature.table,down.genes.id,all.x=T)
		
				down.genes.temp<-cbind(probeid.down.markers,0,1)
			}
		
			}   # marker selection: t-test END
		
			# combine up & down (poor & good) templates or make templates using statistics
		
			if (marker.selection=="ttest" | marker.selection=="snr"){
			pred.template<-NULL
			if (is.null(up.genes.temp)&is.null(down.genes.temp)){
				pred.template<-NULL
			}else{
				pred.template<-as.data.frame(rbind(up.genes.temp,down.genes.temp))
				colnames(pred.template)<-c("probeid","up","down")
			}
			}
			if (marker.selection=="cox.score" | marker.selection=="samr.surv"){
			pred.template<-NULL
			if (is.null(poor.genes.temp)&is.null(good.genes.temp)){
				pred.template<-NULL
			}else{
				pred.template<-as.data.frame(rbind(poor.genes.temp,good.genes.temp))
				colnames(pred.template)<-c("probeid","poor","good")
			}
			}
		
			if (temp.nn.wt=="T"){
			pred.template[,2] <- marker.stat
			pred.template[,3] <- - marker.stat
			}
		#print(pred.template[1:3,])
		
			# prediction for LOO sample
			if (length(pred.template)!=0){
		#print(loo.sample.data[1:3,])
			featured.loo.sample.data<-merge(pred.template,loo.sample.data,sort=F)
			# col 1,2,3,4: probeid,clsUP,clsDOWN,row.normed.data
		
			if (marker.selection=="ttest" | marker.selection=="snr"){
				up.template <- as.numeric(featured.loo.sample.data[,2])
				down.template <- as.numeric(featured.loo.sample.data[,3])
				sample.data <- as.numeric(as.vector(featured.loo.sample.data[,4]))
		#print(paste("up temp: ",up.template,sep=""))
		#print(paste("down temp: ",down.template,sep=""))
		#print(paste("sample data: ",sample.data,sep=""))
		#print(paste("marker stat: ",marker.stat,sep=""))
		
				if (temp.nn.wt=="T"){           # weight by gene ranking statistic
				sample.data <- sample.data*abs(marker.stat)
				}
		#print(paste("wt.sample data: ",sample.data,sep=""))
		
				if (dist.selection=="cosine"){
				distance.to.up<-cosine(up.template,sample.data)
				distance.to.down<-cosine(down.template,sample.data)
				}
				if (dist.selection=="correlation"){
				distance.to.up<-cor(up.template,sample.data)
				distance.to.down<-cor(down.template,sample.data)
				}
		
		#        distance.to.up<-cosine(as.numeric(featured.loo.sample.data[,2]),as.numeric(featured.loo.sample.data[,4]))
		#        distance.to.down<-cosine(as.numeric(featured.loo.sample.data[,3]),as.numeric(featured.loo.sample.data[,4]))
		
				if (distance.to.up>distance.to.down){
				pred.result[loo,3]<-cls.label.1
				pred.result[loo,4]<-distance.to.up
				}
				if (distance.to.up<distance.to.down){
				pred.result[loo,3]<-cls.label.2
				pred.result[loo,4]<-distance.to.down
				}
				if (distance.to.up==distance.to.down){
				pred.result[loo,3]<-"border"
				pred.result[loo,4]<-distance.to.up
				}
			}
		
			if (marker.selection=="cox.score" | marker.selection=="samr.surv"){
				poor.template <- as.numeric(featured.loo.sample.data[,2])
				good.template <- as.numeric(featured.loo.sample.data[,3])
				sample.data <- as.numeric(as.vector(featured.loo.sample.data[,4]))
		#print(paste("poor temp: ",poor.template,sep=""))
		#print(paste("good temp: ",good.template,sep=""))
		#print(paste("sample data: ",sample.data,sep=""))
		#print(paste("marker stat: ",marker.stat,sep=""))
		
				if (temp.nn.wt=="T"){           # weight by gene ranking statistic
				sample.data <- sample.data*abs(marker.stat)
				}
		#print(paste("wt.sample data: ",sample.data,sep=""))
		
				if (dist.selection=="cosine"){
				distance.to.poor<-cosine(poor.template,sample.data)
				distance.to.good<-cosine(good.template,sample.data)
				}
				if (dist.selection=="correlation"){
				distance.to.poor<-cor(poor.template,sample.data)
				distance.to.good<-cor(good.template,sample.data)
				}
		#print(distance.to.poor)
		#print(distance.to.good)
				if (distance.to.poor>distance.to.good){
				pred.result[loo,3]<-"poor"
				pred.result[loo,4]<-distance.to.poor
				}
				if (distance.to.poor<distance.to.good){
				pred.result[loo,3]<-"good"
				pred.result[loo,4]<-distance.to.good
				}
				if (distance.to.poor==distance.to.good){
				pred.result[loo,3]<-"border"
				pred.result[loo,4]<-distance.to.poor
				}
			}
		
			}else{
			pred.result[loo,3]<-"NoMarker"
			pred.result[loo,4]<-NA
			}
		
			# null distribution
		
			if (marker.selection=="cox.score" | marker.selection=="samr.surv"){
			for ( p in 1:nresample){
				if (within.sig=="T"){
				resample.loo.data <- sample(sample.data,length(pred.template[,1]))
				}else{
				resample.loo.data <- as.numeric(sample(loo.sample.data[,2],length(pred.template[,1])))
				}
		#print(resample.loo.data[1:3])
				if (dist.selection=="cosine"){
				resample.stat[loo,p]<-cosine(poor.template,resample.loo.data)
				resample.stat[loo,(nresample+p)]<-cosine(good.template,resample.loo.data)
				}
				if (dist.selection=="correlation"){
				resample.stat[loo,p]<-cor(poor.template,resample.loo.data)
				resample.stat[loo,(nresample+p)]<-cor(good.template,resample.loo.data)
				}
			}
			}
		
			if (marker.selection=="ttest" | marker.selection=="snr"){
			for ( p in 1:nresample){
				if (within.sig=="T"){
				resample.loo.data <- sample(sample.data,length(pred.template[,1]))
				}else{
				resample.loo.data <- as.numeric(sample(loo.sample.data[,2],length(pred.template[,1])))
				}
		#print(resample.loo.data[1:3])
				if (dist.selection=="cosine"){
				resample.stat[loo,p]<-cosine(up.template,resample.loo.data)
				resample.stat[loo,(nresample+p)]<-cosine(down.template,resample.loo.data)
				}
				if (dist.selection=="correlation"){
				resample.stat[loo,p]<-cor(up.template,resample.loo.data)
				resample.stat[loo,(nresample+p)]<-cor(down.template,resample.loo.data)
				}
			}
			}
		
		#print(resample.stat[loo,1:5])
		
			if (histgram.cosine=="T"){
			png(paste("histgram.cosine_",sample.names[loo],".png",sep=""))
				hist(resample.stat[loo,])
			dev.off()
			}
		}  # end LOOCV loop
		
		# confidence of prediction
		orig.resample.stat <- as.data.frame(cbind(as.numeric(as.vector(pred.result[,4])),resample.stat))
		#print(dim(orig.resample.stat))
		rank.orig.resample.stat <- apply(orig.resample.stat,1,rank)
		#print(rank.orig.resample.stat[1,])
		nominal.p <- ((nresample*2)+1-as.vector(rank.orig.resample.stat[1,]))/(nresample*2)
		#  nominal.p <- (nresample+1-as.vector(rank.orig.resample.stat[1,]))/nresample
		#print(nominal.p)
		rank.nominal.p <- rank(nominal.p)
		BH.FDR <- nominal.p*num.sample/rank.nominal.p
		Bonferroni <- nominal.p*num.sample
		BH.FDR[BH.FDR>1] <- 1
		Bonferroni[Bonferroni>1] <- 1
		
		pred.result[,5] <- nominal.p
		pred.result[,6] <- BH.FDR
		pred.result[,7] <- Bonferroni
		
		# count feature usage
		
		if (marker.selection=="ttest" | marker.selection=="snr"){
			count.up<-vector(length=num.gene,mode="numeric")
			count.down<-vector(length=num.gene,mode="numeric")
			percent.up<-vector(length=num.gene,mode="numeric")
			percent.down<-vector(length=num.gene,mode="numeric")
		
			up.feature.matrix<-as.matrix(up.genes.feature.table[,3:length(up.genes.feature.table[1,])])
			up.feature.matrix[up.feature.matrix=="up"]<-1
			down.feature.matrix<-as.matrix(down.genes.feature.table[,3:length(down.genes.feature.table[1,])])
			down.feature.matrix[down.feature.matrix=="down"]<-1
		
			for (i in 1:num.gene){
			count.up[i]<-sum(as.numeric(up.feature.matrix[i,]),na.rm=T)
			percent.up[i]<-count.up[i]/num.sample
			count.down[i]<-sum(as.numeric(down.feature.matrix[i,]),na.rm=T)
			percent.down[i]<-count.down[i]/num.sample
			}
		#    each.down<-down.genes.feature.table[i,3:length(down.genes.feature.table[1,])]
		#    count.down[i]<-length(each.down[each.down=="down"])
		#    percent.down[i]<-count.down[i]/num.sample
		
		
			up.genes.feature.table<-cbind(up.genes.feature.table,count.up,percent.up)
			down.genes.feature.table<-cbind(down.genes.feature.table,count.down,percent.down)
		}
		
		if (marker.selection=="cox.score" | marker.selection=="samr.surv"){
			count.poor<-vector(length=num.gene,mode="numeric")
			count.good<-vector(length=num.gene,mode="numeric")
			percent.poor<-vector(length=num.gene,mode="numeric")
			percent.good<-vector(length=num.gene,mode="numeric")
		
			poor.feature.matrix<-as.matrix(poor.genes.feature.table[,3:length(poor.genes.feature.table[1,])])
			poor.feature.matrix[poor.feature.matrix=="poor"]<-1
			good.feature.matrix<-as.matrix(good.genes.feature.table[,3:length(good.genes.feature.table[1,])])
			good.feature.matrix[good.feature.matrix=="good"]<-1
		
			for (i in 1:num.gene){
			count.poor[i]<-sum(as.numeric(poor.feature.matrix[i,]),na.rm=T)
			percent.poor[i]<-count.poor[i]/num.sample
			count.good[i]<-sum(as.numeric(good.feature.matrix[i,]),na.rm=T)
			percent.good[i]<-count.good[i]/num.sample
			}
		#    each.good<-good.genes.feature.table[i,3:length(good.genes.feature.table[1,])]
		#    count.good[i]<-length(each.good[each.good=="good"])
		#    percent.good[i]<-count.good[i]/num.sample
		
		
			poor.genes.feature.table<-cbind(poor.genes.feature.table,count.poor,percent.poor)
			good.genes.feature.table<-cbind(good.genes.feature.table,count.good,percent.good)
		}
		
		# output
		
		if (marker.selection=="ttest" | marker.selection=="snr"){
			up.features.output.name <- paste(output.file.2,"_",cls.label.1,".features.txt",sep="")
			down.features.output.name <- paste(output.file.2,"_",cls.label.2,".features.txt",sep="")
		
			pred.result.output.name <- paste(output.file.2,"_prediction.result.txt",sep="")
		
			write.table(up.genes.feature.table,up.features.output.name,quote=F,sep="\t",row.names=F)
			write.table(down.genes.feature.table,down.features.output.name,quote=F,sep="\t",row.names=F)
		
			write.table(pred.result,pred.result.output.name,quote=F,sep="\t",row.names=F)
		}
		if (marker.selection=="cox.score" | marker.selection=="samr.surv"){
			poor.features.output.name <- paste(output.file.2,"_poor.features.txt",sep="")
			good.features.output.name <- paste(output.file.2,"_good.features.txt",sep="")
		
			pred.result.output.name <- paste(output.file.2,"_prediction.result.txt",sep="")
		
			write.table(poor.genes.feature.table,poor.features.output.name,quote=F,sep="\t",row.names=F)
			write.table(good.genes.feature.table,good.features.output.name,quote=F,sep="\t",row.names=F)
		
			write.table(pred.result,pred.result.output.name,quote=F,sep="\t",row.names=F)
		}
		
		# end of main
	}	
	
	###########################################################
	#
	#
	# Process results of SurvivalGene_GP.R and loocv.nn_cox_GP.R
	#
	# 
	############################################################
	library(data.table)
	
	df.cox <- fread(file=paste0(output.file,".txt"), header=T, select = c("probeid","statistic"))
	df.good <- fread(file=paste0(output.file,"_LOOCV_Survival_NTP_good.features.txt"), header=T, select = c("probeid","percent.good"))
	df.poor <- fread(file=paste0(output.file,"_LOOCV_Survival_NTP_poor.features.txt"), header=T, select = c("probeid","percent.poor"))
	
	df.good.genes <- df.good[which(df.good$percent.good == 1),]
	df.poor.genes <- df.poor[which(df.poor$percent.poor == 1),]
	
	df.result <- rbind(merge(df.poor.genes,df.cox),merge(df.good.genes,df.cox), fill=T)
	df.result <- df.result[,c(1,3)]
	df.result.ordered <- df.result[order(df.result[,2], decreasing=T),]
	
	#write results .txt files
	write.table(df.result, file=paste0(output.file,"_result.txt"), quote=F, row.names=F, sep="\t")
	write.table(df.result.ordered, file=paste0(output.file,"_result_ordered.txt"), quote=F, row.names=F, sep="\t")

}

