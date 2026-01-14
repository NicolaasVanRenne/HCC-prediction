#filter signature for coefficient of variation
#################################################
#set working directory
#setwd("") 

#set input files
input.exp.filename      = "data_files/collapsed_data/GSE237330_training_RPM.gct"    #gct file of gene expression 
input.features.filename = "output/training_RPM_filtered_signature_result_ordered.txt" #.txt file with weighted genes 
#input.features.filename = "output/2_calculate_signature/training_56_RPM_filtered_signature_result_ordered_FILTER.txt" #.txt file with weighted genes 

#calculate coefficient of variation
exp.dataset<-read.delim(input.exp.filename,header=T,skip=2,check.names=F)
colnames(exp.dataset)[1:2] <- c("ProbeID","GeneName")
exp.df <- exp.dataset[,c(3:ncol(exp.dataset))]
rownames(exp.df) <- exp.dataset$ProbeID

cv_gene <- apply(exp.df, 1, function(x) sd(x) / mean(x))

#plot
hist(cv_gene)

#adapt input features to the .txt file necessary for the genepattern algorithm
#marker genes annotations with 4 columns, Probe ID, Gene name, Class (1,2,...), Weight (optional)
features<-read.delim(input.features.filename,header=T,check.names=F)
if(ncol(features)>2){print("error!!!! more than 2 columns. Set adapt.genepattern == FALSE ????")}
features.cls <- c() ; features.cls[which(features[,2] > 1)] <- 1 ;  features.cls[which(features[,2] < 1)] <- 2
features <- cbind(features[,1],features[,1], features.cls, features[,2])
colnames(features) <- c("ID","description","class","weight")
features <- as.data.frame(features)

cv <- cv_gene[names(cv_gene) %in% features$ID]
features.backup <- features


#filter features for weight (Cox Score)
cox.limit = 0.00 # set limit of Cox scores (absolute values)

features <- features.backup[abs(as.numeric(features.backup$weight)) > cox.limit ,]
nrow(features)


#plot
hist(cv)


# Reorder cv_selection according to the order of IDs in the features dataframe
ordered_cv_selection <- cv[match(features$ID, names(cv))]

# Add the ordered cv_selection as a new column to the features dataframe
features$cv <- ordered_cv_selection

#write table
#write.table(features, file="output/features_cv.txt", row.names=F, quote=F, sep="\t")

#make a dot plot of cox score (weight) and cv
library(ggplot2)
ggplot(features, aes(x = as.numeric(weight), y = as.numeric(cv))) +
  geom_point() +  # Add points for each observation
  theme_minimal() 

#make a bar graph for publication
highlight_data <- features[features$ID %in% c("IGHA1", "IGHA2"), ]
bin_width = 0.25
p=ggplot(features, aes(x = cv)) +
  geom_histogram(binwidth = bin_width, boundary = 0, fill="lightgrey", color="black") +
  #scale_x_continuous(breaks = seq(0, max(features$cv), by = bin_width)) +
  scale_x_continuous(breaks = seq(floor(min(features$cv)), ceiling(max(features$cv)), by = 1)) + # Labels only at integer values
  scale_y_continuous(expand = c(0, 0)) + # This will make the bars touch the Y-axis
  labs(x = "Coefficient of variation", y = "Number of genes", title = "Distribution of Coefficient of variation") +
  theme_classic() +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "blue", size = 1.5)  # Adds the dashed line
#annotate("text", x = max(cv_df$cv) - 1, y = 45, 
#        label = "Reduced signature", size = 6, fontface = "bold") +
#annotate("text", x = max(cv_df$cv) - 1, y = 35, 
#        label = "61 high risk genes\n 93 low risk genes", size = 6)
#geom_text(data = highlight_data, aes(x = cv, y = after_stat(count), 
#         label = paste("cv =", rownames(highlight_data), "(", round(cv,digits=2), ")", sep = "")), 
#        stat = "bin", vjust = 0.5, color = "red", size = 5)
p
#tiff("reduced_hist_cv0.75.tiff", width = 350, height = 350)
#p
#dev.off()


#set limit for cv
cv.limit = 1.25 # set limit of coefficient of variation

features <- features[features$cv > cv.limit ,]
features <- features[ , -5]
features.temp <- features

table(features$class)
sum(table(features$class))

#turn it into gene sets
genesets <- list()
genesets[[1]] <- features[,"ID"][features[,"weight"] > 0] ; names(genesets)[1] <- "High-risk genes"
genesets[[2]] <- features[,"ID"][features[,"weight"] < 0] ; names(genesets)[2] <- "Low-risk genes"

# save
#write.table(features[,c("ID","weight")], file="reduced_signature.txt", sep="\t", row.names = F, col.names = F)
orig_sign <- read.table(file="data_files/gmt_files/signatures.gmt", sep="\t", fill=T)
new_sign <- orig_sign[c(3,4),]
poor.genes <- features[features$weight > 0, 1]
good.genes <- features[features$weight < 0, 1]
new_sign[3,] <- c(rep("Reduced_VanRenne_poor",2), genesets[[1]], rep("", ncol(orig_sign)-2-length(genesets[[1]])))
new_sign[4,] <- c(rep("Reduced_VanRenne_good",2), genesets[[2]], rep("", ncol(orig_sign)-2-length(genesets[[2]])))
#write.table(new_sign, file="data_files/gmt_files/reduced_signature.gmt", sep="\t", row.names = F, col.names = F)



######################################################
# 
# Nearest Template prediction 
# 
#######################################################
#Note: The NTP algorithm was created by Yujin Hoshida. 
#Note: If you use this algorithm, Please cite Hoshida et al. PLoS One 2010, PMID:21124904.
#Note: The algorithm and documentation are available on genepattern (https://www.genepattern.org/)
#Note: Since these algorithms are redistributed from GenePattern, the original copyright notice is included here:

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




#load libraries
library(survival)
library(survminer)
library(patchwork)
library(dplyr)

save.plots=FALSE

#set working directory
#setwd("") 


#set random seed version to old (R version < 3.6) or new (R version > 3.6)
choose.seed = "old" #"old" or "new": set random seed version to old (R version < 3.6; reproduces GenePattern) or new (R version > 3.6)
if(choose.seed == "old"){RNGkind(sample.kind = "Rounding")}else{} #necessary for reproducing Hoshida NEJM 2008 because set.seed changed since R3.6


### 1. Training set
######################################

#set input files
input.exp.filename      = "data_files/collapsed_data/GSE237330_training_RPM.gct"    #gct file of gene expression 
###del###	input.features.filename = "output/2_calculate_signature/training_56_RPM_filtered_signature_result_ordered.txt" #.txt file with weighted genes 
output.name             = "output/3_NTP/training_NTP"
input.meta.data     	= "data_files/meta_data/training_clinical_HCC.txt"


###del###  #adapt input features to the .txt file necessary for the genepattern algorithm
###del###  	#marker genes annotations with 4 columns, Probe ID, Gene name, Class (1,2,...), Weight (optional)
###del###  	features<-read.delim(input.features.filename,header=T,check.names=F)
###del###  	if(ncol(features)>2){print("error!!!! more than 2 columns. Set adapt.genepattern == FALSE ????")}
###del###  	features.cls <- c() ; features.cls[which(features[,2] > 1)] <- 1 ;  features.cls[which(features[,2] < 1)] <- 2
###del###  	features <- cbind(features[,1],features[,1], features.cls, features[,2])
###del###  	colnames(features) <- c("ID","description","class","weight")
###del###  	features <- as.data.frame(features)


########################################################################
# NTPez.R                                    Jan.27, 2009
#                                     Fixed: 06/18/2012
#   Nearest Template Prediction (NTP) for >=2 classes
#   Yujin Hoshida (Broad Institute)
#
#    input: (1) Gene expression dataset (.gct)
#           (2) Marker genes (.txt)
#                  1st column: Probe ID
#                  2nd column: Gene name
#                  3rd column: Class (1,2,...)
#                  4th column (optional): Weight vector
#                           (score of differential expression, etc)
#    output: Prediction result, Marker genes (.xls), heatmap (.png)
#            Files for GenePattern (.gct, .cls)
########################################################################

{
  # temp.nn distance & row normalize
  dist.selection="cosine" # "correlation" or "cosine"
  temp.nn.wt="T"   # only for 2 cls
  
  # resampling to generate null dist
  nresmpl=1000
  
  # outputs
  GenePattern.output="T"
  
  # Seed
  rnd.seed=7392854
  
  #  suppressWarnings()
  
  # Advanced setting
  
  norm.method="row.std" # "row.std.ref","ratio.ref"
  ref.sample.file=NULL
  within.sig="F"
  plot.nominal.p="F"
  plot.distance="F"
  dchip.output="F"
  signature.heatmap="T" # NVR: default "T", set to "F" to avoid heatmap generation 
  FDR.sample.bar=0.05 # NA if not needed
  plot.FDR="T"           # NVR: default "T", set to "F" to avoid FDR plot generation 
  col.range=3         # SD in heatmap
  heatmap.legend=signature.heatmap
  histgram.null.dist="F" # histgram of null dist for the distance
  hist.br=30
  
  # for GenePattern
  
  nresmpl<-as.numeric(nresmpl)
  col.range<-as.numeric(col.range)
  hist.br<-as.numeric(hist.br)  # bin number for resampled dist histgram
  rnd.seed <- as.numeric(rnd.seed)
  if (FDR.sample.bar!="NA"){
    FDR.sample.bar <- as.numeric(FDR.sample.bar)
    if (is.numeric(FDR.sample.bar)==F){
      stop("### Provide numerical value (0~1) for FDR.sample.bar! ###")
    }
  }
  if (is.null(ref.sample.file)){
  }else{
    if (ref.sample.file=="NULL"){
      ref.sample.file <- NULL
    }
  }
  
  # set random seed
  set.seed(rnd.seed)
  
  ### input ###
  
  # selected features used for prediction
  
  ### own features data frame made first!!!### features<-read.delim(input.features.filename,header=T,check.names=F)
  
  ## file format check
  if (length(features[1,])!=3 & length(features[1,])!=4){
    stop("### Please use features file format! ###")
  }
  if (length(features[1,])<4 & temp.nn.wt=="T"){
    temp.nn.wt <- "F"
  }
  third.col<-rownames(table(features[,3]))
  if (is.na(as.numeric(third.col[1]))){
    stop("### The 3rd column of feature file should be numerical! ###")
  }
  
  feat.col.names<-colnames(features)
  feat.col.names[1:2]<-c("ProbeID","GeneName")
  colnames(features)<-feat.col.names
  
  num.features<-length(features[,1])
  num.cls<-length(table(features[,3]))
  feature.col.num <- length(features[1,])
  
  ord<-seq(1:num.features)
  features<-cbind(ord,features)  # add order column to "features"
  
  # expression data
  ## file format check
  if (regexpr(".gct$",input.exp.filename)==-1){
    stop("### Gene expression data should be .gct format! ###")
  }
  exp.dataset<-read.delim(input.exp.filename,header=T,skip=2,check.names=F)
  colnames(exp.dataset)[1:2] <- c("ProbeID","GeneName")
  
  ## Other dataset's mean & SD for row normalization (optional)
  if (!is.null(ref.sample.file)){
    ref.sample <- read.delim(ref.sample.file,header=T)
    if (dim(ref.sample)[2]!=4 & is.numeric(ref.sample[1,3]) & is.numeric(ref.sample[1,4])){
      stop("### mean & SD file format incorrect! ###")
    }
    colnames(ref.sample)[1:4] <- c("ProbeID","SomeName","mean","sd")
    merged.dataset <- merge(ref.sample,exp.dataset,sort=F)
    
    ref.sample <- merged.dataset[,1:4]
    exp.dataset <- merged.dataset[,c(1,5:dim(merged.dataset)[2])]
  }
  
  ProbeID<-exp.dataset[,1]
  gene.names<-exp.dataset[,2]
  num.samples<-(length(exp.dataset[1,])-2)
  exp.dataset<-exp.dataset[-c(1:2)]
  
  exp.for.sample.names<-read.delim(input.exp.filename,header=F,skip=2)  # read sample names
  sample.names<-as.vector(as.matrix(exp.for.sample.names[1,3:length(exp.for.sample.names[1,])]))
  
  # row normalize
  
  normed.exp.dataset<-exp.dataset
  
  if (norm.method=="row.std"){
    exp.mean <- apply(exp.dataset,1,mean,na.rm=T)
    exp.sd <- apply(exp.dataset,1,sd,na.rm=T)
    normed.exp.dataset<-(exp.dataset-exp.mean)/exp.sd
  }
  if (norm.method=="row.std.ref"){
    if (is.null(ref.sample)){
      stop("### Provide reference sample data! ###")
    }
    exp.mean <- as.numeric(as.vector(ref.sample$mean))
    exp.sd <- as.numeric(as.vector(ref.sample$sd))
    normed.exp.dataset<-(exp.dataset-exp.mean)/exp.sd
  }
  if (norm.method=="ratio.ref"){
    if (is.null(ref.sample)){
      stop("### Provide reference sample data! ###")
    }
    exp.mean <- as.numeric(as.vector(ref.sample$mean))
    normed.exp.dataset<- exp.dataset/exp.mean
  }
  
  normed.exp.dataset<-cbind(ProbeID,normed.exp.dataset)
  
  # extract features from normed.exp.dataset
  
  exp.dataset.extract<-merge(features,normed.exp.dataset,sort=F)
  if (length(exp.dataset.extract[,1])<1){
    stop("### No matched probes! ###")
  }
  
  order.extract<-order(exp.dataset.extract[,2])
  exp.dataset.extract<-exp.dataset.extract[order.extract,]
  order.extract.after<-exp.dataset.extract[,2]
  exp.dataset.extract<-exp.dataset.extract[-2]
  
  if (temp.nn.wt=="F"){
    features.extract<-exp.dataset.extract[,1:3]
    if (feature.col.num==4){
      exp.dataset.extract <- exp.dataset.extract[-4]
    }
    features.extract<-cbind(order.extract.after,features.extract) # order:ProbeID:gene name:cls:wt(if any)
    num.features.extract<-length(features.extract[,1])
    
    ProbeID.extract<-as.vector(exp.dataset.extract[,1])
    exp.dataset.extract<-exp.dataset.extract[-c(1:3)]
    rownames(exp.dataset.extract)<-ProbeID.extract
  }
  
  #  temp.nn.wt.vector <- rep(1,num.features)
  
  if (temp.nn.wt=="T" & num.cls==2){
    features.extract<-exp.dataset.extract[,1:4]
    features.extract<-cbind(order.extract.after,features.extract) # order:ProbeID:gene name:cls:wt(if any)
    
    #    if (is.numeric(features[,4])){
    temp.nn.wt.vector <- as.numeric(as.vector(features.extract[,5]))
    #    }else{
    if (is.numeric(temp.nn.wt.vector)==F){
      stop("# Please use numeric values in 4th column!#")
    }
    
    num.features.extract<-length(features.extract[,1])
    
    ProbeID.extract<-as.vector(exp.dataset.extract[,1])
    exp.dataset.extract<-exp.dataset.extract[-c(1:4)]
    rownames(exp.dataset.extract)<-ProbeID.extract
  }
  
  # make template
  
  for (i in 1:num.cls){
    temp.temp<-as.numeric(as.vector(features.extract[,4]))
    temp.temp[temp.temp!=i]<-0
    temp.temp[temp.temp==i]<-1
    eval(parse(text=paste("temp.",i,"<-temp.temp",sep="")))
    #    eval(parse(text=paste("temp\.",i,"<-temp\.temp",sep="")))  ### for < R-2.4.0
  }
  
  # weighted template (only for 2cls)
  
  if (temp.nn.wt=="T" & num.cls==2){
    temp.1 <- temp.nn.wt.vector
    temp.2 <- -temp.nn.wt.vector
  }
  
  ### compute distance and p-value ###
  
  predict.label<-vector(length=num.samples,mode="numeric")
  dist.to.template<-vector(length=num.samples,mode="numeric")
  dist.to.cls1<-vector(length=num.samples,mode="numeric")
  
  rnd.feature.matrix<-matrix(0,nrow=num.features.extract,ncol=nresmpl)
  
  perm.dist.vector<-vector(length=nresmpl*num.cls,mode="numeric")
  nominal.p<-vector(length=num.samples,mode="numeric")
  BH.FDR<-vector(length=num.samples,mode="numeric")
  Bonferroni.p<-vector(length=num.samples,mode="numeric")
  
  for (i in 1:num.samples){
    
    print(paste("sample # ",i,sep=""))
    
    current.sample <- as.vector(exp.dataset.extract[,i])
    
    # compute original distance
    
    orig.dist.to.all.temp <- vector(length=num.cls,mode="numeric")
    
    if (temp.nn.wt=="T"){   # weight sample data
      current.sample <- current.sample*abs(temp.nn.wt.vector)
    }
    
    if (dist.selection=="cosine"){
      for (o in 1:num.cls){      # compute distance to all templates
        eval(parse(text=paste("current.temp <- temp.",o,sep="")))
        #        eval(parse(text=paste("current\.temp <- temp\.",o,sep="")))  ### for < R-2.4.0
        orig.dist.to.all.temp[o]<-sum(current.temp*current.sample)/
          (sqrt(sum(current.temp^2))*sqrt(sum(current.sample^2)))
      }
    }
    if (dist.selection=="correlation"){
      for (o in 1:num.cls){      # compute distance to all templates
        eval(parse(text=paste("current.temp <- temp.",o,sep="")))
        #        eval(parse(text=paste("current\.temp <- temp\.",o,sep="")))  ### for < R-2.4.0
        orig.dist.to.all.temp[o] <- cor(current.temp,current.sample,method="pearson",use="complete.obs")
      }
    }
    
    if (num.cls==2){           # find nearest neighbor (2 classes)
      if (orig.dist.to.all.temp[1]>=orig.dist.to.all.temp[2]){
        predict.label[i]<-1
        dist.to.template[i]<-1-orig.dist.to.all.temp[1]
        dist.to.cls1[i]<--(orig.dist.to.all.temp[1]+1)
      }
      if (orig.dist.to.all.temp[1]<orig.dist.to.all.temp[2]){
        predict.label[i]<-2
        dist.to.template[i]<-1-orig.dist.to.all.temp[2]
        dist.to.cls1[i]<-orig.dist.to.all.temp[2]+1
      }
    }
    
    if (num.cls>2){
      for (o in 1:num.cls){       # find nearest neighbor (>2 classes)
        if (is.na(orig.dist.to.all.temp[o])!=T){
          if (orig.dist.to.all.temp[o]==max(orig.dist.to.all.temp,na.rm=T)){
            predict.label[i]<-o
            dist.to.template[i]<-1-orig.dist.to.all.temp[o]
            dist.to.cls1[i]<-(1-orig.dist.to.all.temp[o])+o
          }
        }
      }
    }
    
    # permutation test
    
    if (within.sig=="F"){     # generate resampled features from all probes
      for (p in 1:nresmpl){
        rnd.feature.matrix[,p]<-sample(normed.exp.dataset[,(i+1)],num.features.extract,replace=F)
      }
    }
    if (within.sig=="T"){     # generate resampled features from only signature genes
      for (p in 1:nresmpl){
        rnd.feature.matrix[,p]<-sample(exp.dataset.extract[,i],num.features.extract,replace=F)
      }
    }
    
    if (temp.nn.wt=="T" & num.cls==2){
      rnd.feature.matrix <- rnd.feature.matrix*abs(temp.nn.wt.vector)
    }
    
    # compute distance to all templates
    if (dist.selection=="cosine"){          # cosine
      for (res in 1:num.cls){
        eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))
        
        prod.sum<-apply(t(t(rnd.feature.matrix)*temp.resmpl),2,sum)
        
        data.sq.sum<-apply(rnd.feature.matrix^2,2,sum)
        temp.sq.sum<-sum(temp.resmpl^2)
        
        perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)]<-
          (1-(prod.sum/(sqrt(data.sq.sum)*sqrt(temp.sq.sum))))
      }
    }
    
    if (dist.selection=="correlation"){          # correlation
      for (res in 1:num.cls){
        eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))
        perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)]<-
          (1-as.vector(cor(rnd.feature.matrix,temp.resmpl,method="pearson",use="complete.obs")))
      }
    }
    
    # compute nominal p-value
    
    combined.stats.rank<-rank(c(dist.to.template[i],perm.dist.vector))
    nominal.p[i]<-combined.stats.rank[1]/length(combined.stats.rank)
    
    # histogram of combined null distributions
    
    if (histgram.null.dist=="T" & capabilities("png")==T){
      png(paste("resampled_",dist.selection,"_dist_histgram_",sample.names[i],".png",sep=""))
      hist(c(dist.to.template[i],perm.dist.vector),br=hist.br,main=paste(sample.names[i],", # resampling: ",nresmpl,sep=""))
      dev.off()
    }
    
  } # main sample loop END
  
  # MCT correction
  
  BH.FDR<-nominal.p*num.samples/rank(nominal.p)
  Bonferroni.p<-nominal.p*num.samples
  
  BH.FDR[BH.FDR>1]<-1
  Bonferroni.p[Bonferroni.p>1]<-1
  
  ### output ###
  
  # prediction results
  
  dist.to.cls1.rank <- rank(dist.to.cls1)
  pred.summary <- cbind(sample.names,predict.label,dist.to.template,dist.to.cls1.rank,
                        nominal.p,BH.FDR,Bonferroni.p)
  
  write.table(pred.summary,paste(output.name,"_prediction_result.txt",sep=""),
              quote=F,sep="\t",row.names=F)
  
  # extracted features
  
  if (temp.nn.wt=="T" & num.cls==2){
    write.table(features.extract[,2:5],paste(output.name,"_features.txt",sep=""),
                quote=F,sep="\t",row.names=F)
  }
  if (temp.nn.wt=="F"){
    write.table(features.extract[,2:4],paste(output.name,"_features.txt",sep=""),
                quote=F,sep="\t",row.names=F)
  }
  
  # sorted exp dataset for heatmap (row normalized)
  
  t.dataset<-t(exp.dataset.extract)               # sort samples
  t.dataset<-cbind(dist.to.cls1,t.dataset)
  ts.dataset<-t.dataset[order(t.dataset[,1]),]
  to.dataset.out<-ts.dataset[,2:(num.features.extract+1)]
  
  sorted.dataset<-t(to.dataset.out)
  heatmap.dataset<-as.matrix(sorted.dataset)
  #  if (.Platform$OS.type == "windows") {
  #    sorted.dataset<-matrix(gsub(" ","",sorted.dataset),ncol=(num.samples+2))
  #  }
  
  # sorted exp dataset for spreadsheets (not normalized)
  
  exp.dataset.gannot<-cbind(ProbeID,exp.dataset)
  exp.dataset.extract<-merge(features.extract,exp.dataset.gannot,sort=F) # redefine exp.dataset.extract
  
  order.extract<-order(exp.dataset.extract[,2])   # sort genes
  exp.dataset.extract<-exp.dataset.extract[order.extract,]
  if (temp.nn.wt=="T" & num.cls==2){
    exp.dataset.extract<-exp.dataset.extract[-c(1:5)]
  }
  if (temp.nn.wt=="F"){
    exp.dataset.extract<-exp.dataset.extract[-c(1:4)]
  }
  
  t.dataset<-t(exp.dataset.extract)               # sort samples
  t.dataset<-cbind(dist.to.cls1,t.dataset)
  ts.dataset<-t.dataset[order(t.dataset[,1]),]
  to.dataset.out<-ts.dataset[,2:(num.features.extract+1)]
  
  sorted.dataset<-t(to.dataset.out)
  sorted.dataset<-cbind(features.extract[,2:3],sorted.dataset)
  
  sorted.dataset.header<-c("ProbeID","GeneName",sample.names[order(t.dataset[,1])])
  sorted.dataset<-t(cbind(sorted.dataset.header,t(sorted.dataset)))
  
  if (.Platform$OS.type == "windows") {
    sorted.dataset<-matrix(gsub(" ","",sorted.dataset),ncol=(num.samples+2))
  }
  
  # output for GenePattern
  
  if (GenePattern.output=="T"){
    
    # exp data
    write.table("#1.2",paste(output.name,"_sorted_dataset.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
    write.table(paste(num.features.extract,num.samples,sep="\t"),paste(output.name,"_sorted_dataset.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table(sorted.dataset,paste(output.name,"_sorted_dataset.gct",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F,append=T)
    
    # cls ffiles (unsorted, sorted)
    cls.out<-matrix(0,nrow=3,ncol=1)
    
    ## unsorted cls
    cls.out[1,]<-paste(num.samples," ",num.cls," 1",sep="")  # line 1
    cls.out[2,]<-paste("# ",paste(unique(predict.label),collapse=" ")) # line 2
    predict.label.out<-as.numeric(predict.label)-1                      # line 3
    cls.out[3,]<-paste(predict.label.out,collapse=" ")
    
    write.table(cls.out,paste(output.name,"_predicted_unsorted.cls",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F)
    
    ## sorted cls
    sorted.predict.label<-sort(as.numeric(predict.label))
    cls.out[2,]<-paste("# ",paste(unique(sorted.predict.label),collapse=" ")) # line 2
    predict.label.out<-sorted.predict.label-1    # line 3
    cls.out[3,]<-paste(predict.label.out,collapse=" ")
    
    write.table(cls.out,paste(output.name,"_predicted_sorted.cls",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F)
    
    # sample info
    sample.info<-cbind(sample.names,predict.label)
    write.table(sample.info,paste(output.name,"_sample_info.txt",sep=""),
                quote=F,sep="\t",row.names=F)
  }
  
  # output for dChip
  
  if (dchip.output=="T"){
    
    # exp data
    write.table(sorted.dataset,paste(output.name,"_dChip_sorted.dataset.txt",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F)
    
    # sample info
    sample.info<-cbind(sample.names,sample.names,predict.label)
    write.table(sample.info,paste(output.name,"_dChip_sample_info.txt",sep=""),
                quote=F,sep="\t",row.names=F)
    
    # gene info
    dchip.gene.info<-cbind(features.extract[,2:3],NA,features.extract[,3],features.extract[,4],NA)
    colnames(dchip.gene.info)<-c("Probe Set","Identifier","FirstOfLocuslink","FirstOfName","FirstOfFunction","Description")
    write.table(dchip.gene.info,paste(output.name,"_dChip_gene_info.txt",sep=""),
                quote=F,sep="\t",row.names=F)
  }
  
  # heatmap
  
  if (signature.heatmap=="T" & capabilities("png")==T){
    
    subclass.col.source <- c("red","blue","yellow","green","purple","orange","lightblue","darkgreen")
    predict.col.vector <- unique(sort(predict.label))
    subclass.col <- subclass.col.source[predict.col.vector]
    
    heatmap.col <- c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000")
    
    heatmap.dataset[heatmap.dataset>col.range] <- col.range
    heatmap.dataset[heatmap.dataset< -col.range] <- -col.range # bug fixed 06/18/2012
    
    ncol.heat <- length(heatmap.dataset[1,])
    nrow.heat <- length(heatmap.dataset[,1])
    
    heatmap.dataset <- apply(heatmap.dataset,2,rev)
    
    num.pred <- as.vector(table(predict.label))
    num.pred.gene <- as.vector(table(features.extract[,4]))
    
    increment.sample <- cumsum(num.pred)
    increment.sample <- c(0,increment.sample)
    
    increment.gene <- cumsum(num.pred.gene)
    increment.gene <- c(0,increment.gene)
    
    png(paste(output.name,"_heatmap.png",sep=""))
    image(1:ncol.heat,1:nrow.heat,t(heatmap.dataset),axes=F,col=heatmap.col,zlim=c(-col.range,col.range),xlim=c(-0.5,(ncol.heat+0.5+round(ncol.heat*0.05))),ylim=c(-0.5,(nrow.heat+0.5+round(nrow.heat*0.08))),xlab=NA,ylab=NA)
    
    for (c in 1:num.cls){                          # gene bar
      rect((ncol.heat+1),0.5,(ncol.heat+0.5+round(ncol.heat*0.05)),(nrow.heat+0.5-increment.gene[c]),col=subclass.col[c],xpd=T,border=F)
    }
    for (c in 1:length(num.pred)){               # sample bar
      rect((0.5+increment.sample[c]),(nrow.heat+2),(ncol.heat+0.5),(nrow.heat+0.5+round(nrow.heat*0.08)),col=subclass.col[c],xpd=T,border=F)
    }
    dev.off()
    
    # heatmap legend
    
    if (heatmap.legend=="T"){
      png(paste(output.name,"_heatmap_legend.png",sep=""))
      par(plt=c(.1,.9,.45,.5))
      a=matrix(seq(1:15),nc=1)
      image(a,col=heatmap.col,xlim=c(0,1),axes=F,yaxt="n")
      box()
      dev.off()
    }
    
    # FDR sample bar
    
    if (FDR.sample.bar!="NA" & num.cls==2){
      fdr.bar.vector <- predict.label[order(dist.to.cls1.rank)]
      fdr.bar.vector[which(BH.FDR[order(dist.to.cls1.rank)]>=FDR.sample.bar)] <- 3
      png(paste(output.name,"_FDR_",FDR.sample.bar,"_sample_bar.png",sep=""))
      par(plt=c(.1,.9,.45,.5))
      a=matrix(fdr.bar.vector,nc=1)
      image(a,col=c("red","blue","gray"),axes=F,yaxt="n")
      dev.off()
    }
    
    if (FDR.sample.bar!="NA" & num.cls>2){
      fdr.bar.vector <- predict.label[order(dist.to.cls1.rank)]
      fdr.bar.vector[which(BH.FDR[order(dist.to.cls1.rank)]>=FDR.sample.bar)] <- (num.cls+1)
      uni.fdr.bar.vector <- sort(unique(fdr.bar.vector))
      if (length(uni.fdr.bar.vector)>1){
        n.sig.cls <- length(uni.fdr.bar.vector)-1
        uni.fdr.bar.vector <- uni.fdr.bar.vector[1:n.sig.cls]
      }else{
        uni.fdr.bar.vector <- NULL
      }
      
      sig.subclass.col <- c(subclass.col.source[1:num.cls],"gray")
      png(paste(output.name,"_FDR_",FDR.sample.bar,"_sample_bar.png",sep=""))
      par(plt=c(.1,.9,.45,.5))
      a=matrix(fdr.bar.vector,nc=1)
      image(a,col=sig.subclass.col,axes=F,yaxt="n")
      dev.off()
    }
    
  }
  
  # plot FDR
  
  if (plot.FDR=="T" & capabilities("png")==T){
    png(paste(output.name,"_FDR.png",sep=""))
    par(plt=c(0.1,0.95,0.4,0.6),las=2)
    plot(BH.FDR[order(dist.to.cls1)],pch=3,col="blue",lwd=2,ylim=c(0,1),main="BH-FDR")
    box(lwd=2)
    dev.off()
  }
  
  # plot nominal-p
  
  if (plot.nominal.p=="T" & capabilities("png")==T){
    png(paste(output.name,"_nominal.p.png",sep=""))
    par(plt=c(0.1,0.95,0.4,0.6),las=2)
    plot(nominal.p[order(dist.to.cls1)],pch=3,col="blue",lwd=2,ylim=c(0,1),main="nominal p-value")
    box(lwd=2)
    dev.off()
  }
  
  # plot distance to template
  
  if (plot.distance=="T" & capabilities("png")==T){
    png(paste(output.name,"_distance.png",sep=""))
    par(plt=c(0.1,0.95,0.4,0.6),las=2)
    plot(dist.to.template[order(dist.to.cls1)],pch=3,col="blue",lwd=2,ylim=c(0,1),main=paste("1 - ",dist.selection,sep=""))
    box(lwd=2)
    dev.off()
  }
  
  
  
  #write .cls file Cls1 (poor) VS Cls2 (good)
  if (GenePattern.output=="T"){
    #sorted
    cls.out.poor.int.good <- cls.out
    cls.out.poor.int.good[1,] <- paste(num.samples," ",max(fdr.bar.vector)," 1",sep="")  # line 1
    if(max(fdr.bar.vector) == 2){cls.out.poor.int.good[2,] <- paste("#", "Poor","Good")}
    if(max(fdr.bar.vector) == 3){cls.out.poor.int.good[2,] <- paste("#", "Poor","Intermediate","Good")}
    cls.out.poor.int.good[3,] <- paste(fdr.bar.vector,collapse=" ")
    write.table(cls.out.poor.int.good,paste(output.name,"_PoorIntGood_predicted_sorted.cls",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
    
    #unsorted
    fdr.bar.vector.unsorted <- predict.label
    fdr.bar.vector.unsorted[which(BH.FDR>=FDR.sample.bar)] <- 3
    cls.out.poor.int.good <- cls.out
    cls.out.poor.int.good[1,] <- paste(num.samples," ",max(fdr.bar.vector.unsorted)," 1",sep="")  # line 1
    if(max(fdr.bar.vector.unsorted) == 2){
      my.cls <- c("#",c("Poor","Good")[unique( fdr.bar.vector.unsorted)])
      cls.out.poor.int.good[2,] <- paste(my.cls, collapse=" ")
    }
    if(max(fdr.bar.vector.unsorted) == 3){
      my.cls <- c("#",c("Poor","Good","Intermediate")[unique( fdr.bar.vector.unsorted)])
      cls.out.poor.int.good[2,] <- paste0(my.cls, collapse=" ")
    }
    cls.out.poor.int.good[3,] <- paste(fdr.bar.vector.unsorted,collapse=" ")
    write.table(cls.out.poor.int.good,paste(output.name,"_PoorIntGood_predicted_unsorted.cls",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
  }
}	

#note: the output will be printed to your working directory
#note: the NTP output itself is stored in the object pred.summary

#how many genes of the signature are in this dataset? 
table(features.extract$class)
t1=table(features.extract$class)

#Stratify samples according to NTP prediction: BH<0.05
prognosis <- rep(0,nrow(pred.summary)) 
for(i in seq_along(prognosis)){
  if(pred.summary[i,2] == 1 & pred.summary[i,6] < 0.05){ prognosis[i] <- "poor"}
  if(pred.summary[i,2] == 2 & pred.summary[i,6] < 0.05){ prognosis[i] <- "good"}
  if(                         pred.summary[i,6] >= 0.05){ prognosis[i] <- "intermediate"}
}

# Create Kaplan Meier curves
#read clinical data		  
meta.data <- read.table(file = input.meta.data, header=T, sep="\t", row.names=1)
km.df <- cbind(select(meta.data, c("status","time")), prognosis)

#set event and time 
set.time  = "time"    #column name of time
set.event = "status"        #column name of event

#create survival objects
poor.int.good        <- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ prognosis, data = km.df)
poor.good.remove.int <- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ prognosis, data = km.df[which(km.df$prognosis != "intermediate"),])

p.pairwise 			<- pairwise_survdiff(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ prognosis, data = km.df, p.adjust.method = "none")

#plot event curve
event.plot.1 <- ggsurvplot(poor.int.good,         fun= "event", pval = surv_pvalue(poor.int.good)$pval,        conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata", linetype = "strata",   palette = c("mediumblue", "darkgrey", "tomato")) #poor VS intermediate VS good 
event.plot.1.no.pvalue <- ggsurvplot(poor.int.good,         fun= "event",        conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata", linetype = "strata",  palette = c("mediumblue", "darkgrey", "tomato")) #poor VS intermediate VS good 
event.plot.2 <- ggsurvplot(poor.good.remove.int,  fun= "event", pval = surv_pvalue(poor.good.remove.int)$pval, conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata", linetype = "strata", palette = c("mediumblue", "tomato")) #poor VS good (remove intermediates)

p1=arrange_ggsurvplots(list(event.plot.1,event.plot.2), ncol=2, title= "Training NTP prediction")		  
p1
#create truncated Kaplan Meier curve when patients at risk drops below 10% for visualization purpose (not for p-value!)
#store p-values
pval.plot.1 <- surv_pvalue(poor.int.good)$pval
pval.plot.2	<- surv_pvalue(poor.good.remove.int)$pval

#settings
cutoff.at.risk.percentage = 10 #set percentage to cut eg. 'cutoff.at.risk = 10' for 10%
#set.event = "status" #set event to find 10% at risk

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
event.plot 				<- ggsurvplot(poor.int.good,  fun= "event", pval = surv_pvalue(poor.int.good)$pval, pval.coord = c(0,1), conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata",  palette = c("mediumblue", "darkgrey", "tomato"), xlim = c(0,cutoff.time.x.axis )) 
event.plot.no.pvalue 	<- ggsurvplot(poor.int.good,  fun= "event",        conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata",  palette = c("mediumblue", "darkgrey", "tomato"), xlim = c(0,cutoff.time.x.axis )) 
event.plot.2 			<- ggsurvplot(poor.good.remove.int,  fun= "event", pval = surv_pvalue(poor.good.remove.int)$pval, pval.coord = c(0,1), conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata", palette = c("mediumblue", "tomato"), xlim = c(0,cutoff.time.x.axis )) 
event.plot.3 			<- ggsurvplot(poor.int.good,  fun= "event", pval = paste0("poor VS good \n p=", surv_pvalue(poor.good.remove.int)$pval), pval.coord = c(0,0.95), conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata", palette = c("mediumblue", "darkgrey","tomato"), xlim = c(0,cutoff.time.x.axis )) 
arrange_ggsurvplots(list(event.plot,event.plot.2), ncol=2, title= "TITLE HERE")		


#publication picture
p=ggsurvplot(poor.int.good,
             fun= "event",
             pval = paste0("poor VS good \np=", signif(surv_pvalue(poor.good.remove.int)$pval,3)),
             pval.coord = c(0,0.95),
             pval.size=7,
             conf.int = FALSE,
             risk.table = TRUE, risk.table.fontsize=6,
             risk.table.height=0.3,
             risk.table.y.text = TRUE,
             font.x=20,
             font.y=20,
             font.tickslab=20,
             palette = c("mediumblue", "darkgrey","tomato"),
             xlim = c(0,cutoff.time.x.axis),
             legend = 'none',
             xlab="Time (years)",
             ylab="HCC incidence"
)
p$table <- p$table + 
  theme(axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) + ylab("Prognosis") + 
  scale_y_discrete(labels=c('poor', 'intermediate', 'good'))
p
p1 <- p
if(save.plots==TRUE){  tiff('output/plots/KM_trainingx.tiff', units="in", width=7, height=7, res=300, compression = 'lzw');print(p);dev.off()  }



### 2. Validation set
######################################

#set input files
input.exp.filename      = "data_files/collapsed_data/GSE237331_validation_RPM.gct"    #gct file of gene expression 
output.name             = "output/3_NTP/validation_NTP"
input.meta.data     	= "data_files/meta_data/validation_clinical_HCC.txt"

#re-load features
features <- features.temp

########################################################################
# NTPez.R                                    Jan.27, 2009
#                                     Fixed: 06/18/2012
#   Nearest Template Prediction (NTP) for >=2 classes
#   Yujin Hoshida (Broad Institute)
#
#    input: (1) Gene expression dataset (.gct)
#           (2) Marker genes (.txt)
#                  1st column: Probe ID
#                  2nd column: Gene name
#                  3rd column: Class (1,2,...)
#                  4th column (optional): Weight vector
#                           (score of differential expression, etc)
#    output: Prediction result, Marker genes (.xls), heatmap (.png)
#            Files for GenePattern (.gct, .cls)
########################################################################

{
  # temp.nn distance & row normalize
  dist.selection="cosine" # "correlation" or "cosine"
  temp.nn.wt="T"   # only for 2 cls
  
  # resampling to generate null dist
  nresmpl=1000
  
  # outputs
  GenePattern.output="T"
  
  # Seed
  rnd.seed=7392854
  
  #  suppressWarnings()
  
  # Advanced setting
  
  norm.method="row.std" # "row.std.ref","ratio.ref"
  ref.sample.file=NULL
  within.sig="F"
  plot.nominal.p="F"
  plot.distance="F"
  dchip.output="F"
  signature.heatmap="T" # NVR: default "T", set to "F" to avoid heatmap generation 
  FDR.sample.bar=0.05 # NA if not needed
  plot.FDR="T"           # NVR: default "T", set to "F" to avoid FDR plot generation 
  col.range=3         # SD in heatmap
  heatmap.legend=signature.heatmap
  histgram.null.dist="F" # histgram of null dist for the distance
  hist.br=30
  
  # for GenePattern
  
  nresmpl<-as.numeric(nresmpl)
  col.range<-as.numeric(col.range)
  hist.br<-as.numeric(hist.br)  # bin number for resampled dist histgram
  rnd.seed <- as.numeric(rnd.seed)
  if (FDR.sample.bar!="NA"){
    FDR.sample.bar <- as.numeric(FDR.sample.bar)
    if (is.numeric(FDR.sample.bar)==F){
      stop("### Provide numerical value (0~1) for FDR.sample.bar! ###")
    }
  }
  if (is.null(ref.sample.file)){
  }else{
    if (ref.sample.file=="NULL"){
      ref.sample.file <- NULL
    }
  }
  
  # set random seed
  set.seed(rnd.seed)
  
  ### input ###
  
  # selected features used for prediction
  
  ### own features data frame made first!!!### features<-read.delim(input.features.filename,header=T,check.names=F)
  
  ## file format check
  if (length(features[1,])!=3 & length(features[1,])!=4){
    stop("### Please use features file format! ###")
  }
  if (length(features[1,])<4 & temp.nn.wt=="T"){
    temp.nn.wt <- "F"
  }
  third.col<-rownames(table(features[,3]))
  if (is.na(as.numeric(third.col[1]))){
    stop("### The 3rd column of feature file should be numerical! ###")
  }
  
  feat.col.names<-colnames(features)
  feat.col.names[1:2]<-c("ProbeID","GeneName")
  colnames(features)<-feat.col.names
  
  num.features<-length(features[,1])
  num.cls<-length(table(features[,3]))
  feature.col.num <- length(features[1,])
  
  ord<-seq(1:num.features)
  features<-cbind(ord,features)  # add order column to "features"
  
  # expression data
  ## file format check
  if (regexpr(".gct$",input.exp.filename)==-1){
    stop("### Gene expression data should be .gct format! ###")
  }
  exp.dataset<-read.delim(input.exp.filename,header=T,skip=2,check.names=F)
  colnames(exp.dataset)[1:2] <- c("ProbeID","GeneName")
  
  ## Other dataset's mean & SD for row normalization (optional)
  if (!is.null(ref.sample.file)){
    ref.sample <- read.delim(ref.sample.file,header=T)
    if (dim(ref.sample)[2]!=4 & is.numeric(ref.sample[1,3]) & is.numeric(ref.sample[1,4])){
      stop("### mean & SD file format incorrect! ###")
    }
    colnames(ref.sample)[1:4] <- c("ProbeID","SomeName","mean","sd")
    merged.dataset <- merge(ref.sample,exp.dataset,sort=F)
    
    ref.sample <- merged.dataset[,1:4]
    exp.dataset <- merged.dataset[,c(1,5:dim(merged.dataset)[2])]
  }
  
  ProbeID<-exp.dataset[,1]
  gene.names<-exp.dataset[,2]
  num.samples<-(length(exp.dataset[1,])-2)
  exp.dataset<-exp.dataset[-c(1:2)]
  
  exp.for.sample.names<-read.delim(input.exp.filename,header=F,skip=2)  # read sample names
  sample.names<-as.vector(as.matrix(exp.for.sample.names[1,3:length(exp.for.sample.names[1,])]))
  
  # row normalize
  
  normed.exp.dataset<-exp.dataset
  
  if (norm.method=="row.std"){
    exp.mean <- apply(exp.dataset,1,mean,na.rm=T)
    exp.sd <- apply(exp.dataset,1,sd,na.rm=T)
    normed.exp.dataset<-(exp.dataset-exp.mean)/exp.sd
  }
  if (norm.method=="row.std.ref"){
    if (is.null(ref.sample)){
      stop("### Provide reference sample data! ###")
    }
    exp.mean <- as.numeric(as.vector(ref.sample$mean))
    exp.sd <- as.numeric(as.vector(ref.sample$sd))
    normed.exp.dataset<-(exp.dataset-exp.mean)/exp.sd
  }
  if (norm.method=="ratio.ref"){
    if (is.null(ref.sample)){
      stop("### Provide reference sample data! ###")
    }
    exp.mean <- as.numeric(as.vector(ref.sample$mean))
    normed.exp.dataset<- exp.dataset/exp.mean
  }
  
  normed.exp.dataset<-cbind(ProbeID,normed.exp.dataset)
  
  # extract features from normed.exp.dataset
  
  exp.dataset.extract<-merge(features,normed.exp.dataset,sort=F)
  if (length(exp.dataset.extract[,1])<1){
    stop("### No matched probes! ###")
  }
  
  order.extract<-order(exp.dataset.extract[,2])
  exp.dataset.extract<-exp.dataset.extract[order.extract,]
  order.extract.after<-exp.dataset.extract[,2]
  exp.dataset.extract<-exp.dataset.extract[-2]
  
  if (temp.nn.wt=="F"){
    features.extract<-exp.dataset.extract[,1:3]
    if (feature.col.num==4){
      exp.dataset.extract <- exp.dataset.extract[-4]
    }
    features.extract<-cbind(order.extract.after,features.extract) # order:ProbeID:gene name:cls:wt(if any)
    num.features.extract<-length(features.extract[,1])
    
    ProbeID.extract<-as.vector(exp.dataset.extract[,1])
    exp.dataset.extract<-exp.dataset.extract[-c(1:3)]
    rownames(exp.dataset.extract)<-ProbeID.extract
  }
  
  #  temp.nn.wt.vector <- rep(1,num.features)
  
  if (temp.nn.wt=="T" & num.cls==2){
    features.extract<-exp.dataset.extract[,1:4]
    features.extract<-cbind(order.extract.after,features.extract) # order:ProbeID:gene name:cls:wt(if any)
    
    #    if (is.numeric(features[,4])){
    temp.nn.wt.vector <- as.numeric(as.vector(features.extract[,5]))
    #    }else{
    if (is.numeric(temp.nn.wt.vector)==F){
      stop("# Please use numeric values in 4th column!#")
    }
    
    num.features.extract<-length(features.extract[,1])
    
    ProbeID.extract<-as.vector(exp.dataset.extract[,1])
    exp.dataset.extract<-exp.dataset.extract[-c(1:4)]
    rownames(exp.dataset.extract)<-ProbeID.extract
  }
  
  # make template
  
  for (i in 1:num.cls){
    temp.temp<-as.numeric(as.vector(features.extract[,4]))
    temp.temp[temp.temp!=i]<-0
    temp.temp[temp.temp==i]<-1
    eval(parse(text=paste("temp.",i,"<-temp.temp",sep="")))
    #    eval(parse(text=paste("temp\.",i,"<-temp\.temp",sep="")))  ### for < R-2.4.0
  }
  
  # weighted template (only for 2cls)
  
  if (temp.nn.wt=="T" & num.cls==2){
    temp.1 <- temp.nn.wt.vector
    temp.2 <- -temp.nn.wt.vector
  }
  
  ### compute distance and p-value ###
  
  predict.label<-vector(length=num.samples,mode="numeric")
  dist.to.template<-vector(length=num.samples,mode="numeric")
  dist.to.cls1<-vector(length=num.samples,mode="numeric")
  
  rnd.feature.matrix<-matrix(0,nrow=num.features.extract,ncol=nresmpl)
  
  perm.dist.vector<-vector(length=nresmpl*num.cls,mode="numeric")
  nominal.p<-vector(length=num.samples,mode="numeric")
  BH.FDR<-vector(length=num.samples,mode="numeric")
  Bonferroni.p<-vector(length=num.samples,mode="numeric")
  
  for (i in 1:num.samples){
    
    print(paste("sample # ",i,sep=""))
    
    current.sample <- as.vector(exp.dataset.extract[,i])
    
    # compute original distance
    
    orig.dist.to.all.temp <- vector(length=num.cls,mode="numeric")
    
    if (temp.nn.wt=="T"){   # weight sample data
      current.sample <- current.sample*abs(temp.nn.wt.vector)
    }
    
    if (dist.selection=="cosine"){
      for (o in 1:num.cls){      # compute distance to all templates
        eval(parse(text=paste("current.temp <- temp.",o,sep="")))
        #        eval(parse(text=paste("current\.temp <- temp\.",o,sep="")))  ### for < R-2.4.0
        orig.dist.to.all.temp[o]<-sum(current.temp*current.sample)/
          (sqrt(sum(current.temp^2))*sqrt(sum(current.sample^2)))
      }
    }
    if (dist.selection=="correlation"){
      for (o in 1:num.cls){      # compute distance to all templates
        eval(parse(text=paste("current.temp <- temp.",o,sep="")))
        #        eval(parse(text=paste("current\.temp <- temp\.",o,sep="")))  ### for < R-2.4.0
        orig.dist.to.all.temp[o] <- cor(current.temp,current.sample,method="pearson",use="complete.obs")
      }
    }
    
    if (num.cls==2){           # find nearest neighbor (2 classes)
      if (orig.dist.to.all.temp[1]>=orig.dist.to.all.temp[2]){
        predict.label[i]<-1
        dist.to.template[i]<-1-orig.dist.to.all.temp[1]
        dist.to.cls1[i]<--(orig.dist.to.all.temp[1]+1)
      }
      if (orig.dist.to.all.temp[1]<orig.dist.to.all.temp[2]){
        predict.label[i]<-2
        dist.to.template[i]<-1-orig.dist.to.all.temp[2]
        dist.to.cls1[i]<-orig.dist.to.all.temp[2]+1
      }
    }
    
    if (num.cls>2){
      for (o in 1:num.cls){       # find nearest neighbor (>2 classes)
        if (is.na(orig.dist.to.all.temp[o])!=T){
          if (orig.dist.to.all.temp[o]==max(orig.dist.to.all.temp,na.rm=T)){
            predict.label[i]<-o
            dist.to.template[i]<-1-orig.dist.to.all.temp[o]
            dist.to.cls1[i]<-(1-orig.dist.to.all.temp[o])+o
          }
        }
      }
    }
    
    # permutation test
    
    if (within.sig=="F"){     # generate resampled features from all probes
      for (p in 1:nresmpl){
        rnd.feature.matrix[,p]<-sample(normed.exp.dataset[,(i+1)],num.features.extract,replace=F)
      }
    }
    if (within.sig=="T"){     # generate resampled features from only signature genes
      for (p in 1:nresmpl){
        rnd.feature.matrix[,p]<-sample(exp.dataset.extract[,i],num.features.extract,replace=F)
      }
    }
    
    if (temp.nn.wt=="T" & num.cls==2){
      rnd.feature.matrix <- rnd.feature.matrix*abs(temp.nn.wt.vector)
    }
    
    # compute distance to all templates
    if (dist.selection=="cosine"){          # cosine
      for (res in 1:num.cls){
        eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))
        
        prod.sum<-apply(t(t(rnd.feature.matrix)*temp.resmpl),2,sum)
        
        data.sq.sum<-apply(rnd.feature.matrix^2,2,sum)
        temp.sq.sum<-sum(temp.resmpl^2)
        
        perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)]<-
          (1-(prod.sum/(sqrt(data.sq.sum)*sqrt(temp.sq.sum))))
      }
    }
    
    if (dist.selection=="correlation"){          # correlation
      for (res in 1:num.cls){
        eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))
        perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)]<-
          (1-as.vector(cor(rnd.feature.matrix,temp.resmpl,method="pearson",use="complete.obs")))
      }
    }
    
    # compute nominal p-value
    
    combined.stats.rank<-rank(c(dist.to.template[i],perm.dist.vector))
    nominal.p[i]<-combined.stats.rank[1]/length(combined.stats.rank)
    
    # histgram of combined null distributions
    
    if (histgram.null.dist=="T" & capabilities("png")==T){
      png(paste("resampled_",dist.selection,"_dist_histgram_",sample.names[i],".png",sep=""))
      hist(c(dist.to.template[i],perm.dist.vector),br=hist.br,main=paste(sample.names[i],", # resampling: ",nresmpl,sep=""))
      dev.off()
    }
    
  } # main sample loop END
  
  # MCT correction
  
  BH.FDR<-nominal.p*num.samples/rank(nominal.p)
  Bonferroni.p<-nominal.p*num.samples
  
  BH.FDR[BH.FDR>1]<-1
  Bonferroni.p[Bonferroni.p>1]<-1
  
  ### output ###
  
  # prediction results
  
  dist.to.cls1.rank <- rank(dist.to.cls1)
  pred.summary <- cbind(sample.names,predict.label,dist.to.template,dist.to.cls1.rank,
                        nominal.p,BH.FDR,Bonferroni.p)
  
  write.table(pred.summary,paste(output.name,"_prediction_result.txt",sep=""),
              quote=F,sep="\t",row.names=F)
  
  # extracted features
  
  if (temp.nn.wt=="T" & num.cls==2){
    write.table(features.extract[,2:5],paste(output.name,"_features.txt",sep=""),
                quote=F,sep="\t",row.names=F)
  }
  if (temp.nn.wt=="F"){
    write.table(features.extract[,2:4],paste(output.name,"_features.txt",sep=""),
                quote=F,sep="\t",row.names=F)
  }
  
  # sorted exp dataset for heatmap (row normalized)
  
  t.dataset<-t(exp.dataset.extract)               # sort samples
  t.dataset<-cbind(dist.to.cls1,t.dataset)
  ts.dataset<-t.dataset[order(t.dataset[,1]),]
  to.dataset.out<-ts.dataset[,2:(num.features.extract+1)]
  
  sorted.dataset<-t(to.dataset.out)
  heatmap.dataset<-as.matrix(sorted.dataset)
  #  if (.Platform$OS.type == "windows") {
  #    sorted.dataset<-matrix(gsub(" ","",sorted.dataset),ncol=(num.samples+2))
  #  }
  
  # sorted exp dataset for spreadsheets (not normalized)
  
  exp.dataset.gannot<-cbind(ProbeID,exp.dataset)
  exp.dataset.extract<-merge(features.extract,exp.dataset.gannot,sort=F) # redefine exp.dataset.extract
  
  order.extract<-order(exp.dataset.extract[,2])   # sort genes
  exp.dataset.extract<-exp.dataset.extract[order.extract,]
  if (temp.nn.wt=="T" & num.cls==2){
    exp.dataset.extract<-exp.dataset.extract[-c(1:5)]
  }
  if (temp.nn.wt=="F"){
    exp.dataset.extract<-exp.dataset.extract[-c(1:4)]
  }
  
  t.dataset<-t(exp.dataset.extract)               # sort samples
  t.dataset<-cbind(dist.to.cls1,t.dataset)
  ts.dataset<-t.dataset[order(t.dataset[,1]),]
  to.dataset.out<-ts.dataset[,2:(num.features.extract+1)]
  
  sorted.dataset<-t(to.dataset.out)
  sorted.dataset<-cbind(features.extract[,2:3],sorted.dataset)
  
  sorted.dataset.header<-c("ProbeID","GeneName",sample.names[order(t.dataset[,1])])
  sorted.dataset<-t(cbind(sorted.dataset.header,t(sorted.dataset)))
  
  if (.Platform$OS.type == "windows") {
    sorted.dataset<-matrix(gsub(" ","",sorted.dataset),ncol=(num.samples+2))
  }
  
  # output for GenePattern
  
  if (GenePattern.output=="T"){
    
    # exp data
    write.table("#1.2",paste(output.name,"_sorted_dataset.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
    write.table(paste(num.features.extract,num.samples,sep="\t"),paste(output.name,"_sorted_dataset.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table(sorted.dataset,paste(output.name,"_sorted_dataset.gct",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F,append=T)
    
    # cls ffiles (unsorted, sorted)
    cls.out<-matrix(0,nrow=3,ncol=1)
    
    ## unsorted cls
    cls.out[1,]<-paste(num.samples," ",num.cls," 1",sep="")  # line 1
    cls.out[2,]<-paste("# ",paste(unique(predict.label),collapse=" ")) # line 2
    predict.label.out<-as.numeric(predict.label)-1                      # line 3
    cls.out[3,]<-paste(predict.label.out,collapse=" ")
    
    write.table(cls.out,paste(output.name,"_predicted_unsorted.cls",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F)
    
    ## sorted cls
    sorted.predict.label<-sort(as.numeric(predict.label))
    cls.out[2,]<-paste("# ",paste(unique(sorted.predict.label),collapse=" ")) # line 2
    predict.label.out<-sorted.predict.label-1    # line 3
    cls.out[3,]<-paste(predict.label.out,collapse=" ")
    
    write.table(cls.out,paste(output.name,"_predicted_sorted.cls",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F)
    
    # sample info
    sample.info<-cbind(sample.names,predict.label)
    write.table(sample.info,paste(output.name,"_sample_info.txt",sep=""),
                quote=F,sep="\t",row.names=F)
  }
  
  # output for dChip
  
  if (dchip.output=="T"){
    
    # exp data
    write.table(sorted.dataset,paste(output.name,"_dChip_sorted.dataset.txt",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F)
    
    # sample info
    sample.info<-cbind(sample.names,sample.names,predict.label)
    write.table(sample.info,paste(output.name,"_dChip_sample_info.txt",sep=""),
                quote=F,sep="\t",row.names=F)
    
    # gene info
    dchip.gene.info<-cbind(features.extract[,2:3],NA,features.extract[,3],features.extract[,4],NA)
    colnames(dchip.gene.info)<-c("Probe Set","Identifier","FirstOfLocuslink","FirstOfName","FirstOfFunction","Description")
    write.table(dchip.gene.info,paste(output.name,"_dChip_gene_info.txt",sep=""),
                quote=F,sep="\t",row.names=F)
  }
  
  # heatmap
  
  if (signature.heatmap=="T" & capabilities("png")==T){
    
    subclass.col.source <- c("red","blue","yellow","green","purple","orange","lightblue","darkgreen")
    predict.col.vector <- unique(sort(predict.label))
    subclass.col <- subclass.col.source[predict.col.vector]
    
    heatmap.col <- c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000")
    
    heatmap.dataset[heatmap.dataset>col.range] <- col.range
    heatmap.dataset[heatmap.dataset< -col.range] <- -col.range # bug fixed 06/18/2012
    
    ncol.heat <- length(heatmap.dataset[1,])
    nrow.heat <- length(heatmap.dataset[,1])
    
    heatmap.dataset <- apply(heatmap.dataset,2,rev)
    
    num.pred <- as.vector(table(predict.label))
    num.pred.gene <- as.vector(table(features.extract[,4]))
    
    increment.sample <- cumsum(num.pred)
    increment.sample <- c(0,increment.sample)
    
    increment.gene <- cumsum(num.pred.gene)
    increment.gene <- c(0,increment.gene)
    
    png(paste(output.name,"_heatmap.png",sep=""))
    image(1:ncol.heat,1:nrow.heat,t(heatmap.dataset),axes=F,col=heatmap.col,zlim=c(-col.range,col.range),xlim=c(-0.5,(ncol.heat+0.5+round(ncol.heat*0.05))),ylim=c(-0.5,(nrow.heat+0.5+round(nrow.heat*0.08))),xlab=NA,ylab=NA)
    
    for (c in 1:num.cls){                          # gene bar
      rect((ncol.heat+1),0.5,(ncol.heat+0.5+round(ncol.heat*0.05)),(nrow.heat+0.5-increment.gene[c]),col=subclass.col[c],xpd=T,border=F)
    }
    for (c in 1:length(num.pred)){               # sample bar
      rect((0.5+increment.sample[c]),(nrow.heat+2),(ncol.heat+0.5),(nrow.heat+0.5+round(nrow.heat*0.08)),col=subclass.col[c],xpd=T,border=F)
    }
    dev.off()
    
    # heatmap legend
    
    if (heatmap.legend=="T"){
      png(paste(output.name,"_heatmap_legend.png",sep=""))
      par(plt=c(.1,.9,.45,.5))
      a=matrix(seq(1:15),nc=1)
      image(a,col=heatmap.col,xlim=c(0,1),axes=F,yaxt="n")
      box()
      dev.off()
    }
    
    # FDR sample bar
    
    if (FDR.sample.bar!="NA" & num.cls==2){
      fdr.bar.vector <- predict.label[order(dist.to.cls1.rank)]
      fdr.bar.vector[which(BH.FDR[order(dist.to.cls1.rank)]>=FDR.sample.bar)] <- 3
      png(paste(output.name,"_FDR_",FDR.sample.bar,"_sample_bar.png",sep=""))
      par(plt=c(.1,.9,.45,.5))
      a=matrix(fdr.bar.vector,nc=1)
      image(a,col=c("red","blue","gray"),axes=F,yaxt="n")
      dev.off()
    }
    
    if (FDR.sample.bar!="NA" & num.cls>2){
      fdr.bar.vector <- predict.label[order(dist.to.cls1.rank)]
      fdr.bar.vector[which(BH.FDR[order(dist.to.cls1.rank)]>=FDR.sample.bar)] <- (num.cls+1)
      uni.fdr.bar.vector <- sort(unique(fdr.bar.vector))
      if (length(uni.fdr.bar.vector)>1){
        n.sig.cls <- length(uni.fdr.bar.vector)-1
        uni.fdr.bar.vector <- uni.fdr.bar.vector[1:n.sig.cls]
      }else{
        uni.fdr.bar.vector <- NULL
      }
      
      sig.subclass.col <- c(subclass.col.source[1:num.cls],"gray")
      png(paste(output.name,"_FDR_",FDR.sample.bar,"_sample_bar.png",sep=""))
      par(plt=c(.1,.9,.45,.5))
      a=matrix(fdr.bar.vector,nc=1)
      image(a,col=sig.subclass.col,axes=F,yaxt="n")
      dev.off()
    }
    
  }
  
  # plot FDR
  
  if (plot.FDR=="T" & capabilities("png")==T){
    png(paste(output.name,"_FDR.png",sep=""))
    par(plt=c(0.1,0.95,0.4,0.6),las=2)
    plot(BH.FDR[order(dist.to.cls1)],pch=3,col="blue",lwd=2,ylim=c(0,1),main="BH-FDR")
    box(lwd=2)
    dev.off()
  }
  
  # plot nominal-p
  
  if (plot.nominal.p=="T" & capabilities("png")==T){
    png(paste(output.name,"_nominal.p.png",sep=""))
    par(plt=c(0.1,0.95,0.4,0.6),las=2)
    plot(nominal.p[order(dist.to.cls1)],pch=3,col="blue",lwd=2,ylim=c(0,1),main="nominal p-value")
    box(lwd=2)
    dev.off()
  }
  
  # plot distance to template
  
  if (plot.distance=="T" & capabilities("png")==T){
    png(paste(output.name,"_distance.png",sep=""))
    par(plt=c(0.1,0.95,0.4,0.6),las=2)
    plot(dist.to.template[order(dist.to.cls1)],pch=3,col="blue",lwd=2,ylim=c(0,1),main=paste("1 - ",dist.selection,sep=""))
    box(lwd=2)
    dev.off()
  }
  
  
  
  #write .cls file Cls1 (poor) VS Cls2 (good)
  if (GenePattern.output=="T"){
    #sorted
    cls.out.poor.int.good <- cls.out
    cls.out.poor.int.good[1,] <- paste(num.samples," ",max(fdr.bar.vector)," 1",sep="")  # line 1
    if(max(fdr.bar.vector) == 2){cls.out.poor.int.good[2,] <- paste("#", "Poor","Good")}
    if(max(fdr.bar.vector) == 3){cls.out.poor.int.good[2,] <- paste("#", "Poor","Intermediate","Good")}
    cls.out.poor.int.good[3,] <- paste(fdr.bar.vector,collapse=" ")
    write.table(cls.out.poor.int.good,paste(output.name,"_PoorIntGood_predicted_sorted.cls",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
    
    #unsorted
    fdr.bar.vector.unsorted <- predict.label
    fdr.bar.vector.unsorted[which(BH.FDR>=FDR.sample.bar)] <- 3
    cls.out.poor.int.good <- cls.out
    cls.out.poor.int.good[1,] <- paste(num.samples," ",max(fdr.bar.vector.unsorted)," 1",sep="")  # line 1
    if(max(fdr.bar.vector.unsorted) == 2){
      my.cls <- c("#",c("Poor","Good")[unique( fdr.bar.vector.unsorted)])
      cls.out.poor.int.good[2,] <- paste(my.cls, collapse=" ")
    }
    if(max(fdr.bar.vector.unsorted) == 3){
      my.cls <- c("#",c("Poor","Good","Intermediate")[unique( fdr.bar.vector.unsorted)])
      cls.out.poor.int.good[2,] <- paste0(my.cls, collapse=" ")
    }
    cls.out.poor.int.good[3,] <- paste(fdr.bar.vector.unsorted,collapse=" ")
    write.table(cls.out.poor.int.good,paste(output.name,"_PoorIntGood_predicted_unsorted.cls",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
  }
}	

#note: the output will be printed to your working directory
#note: the NTP output itself is stored in the object pred.summary

#how many genes of the signature are in this dataset? 
table(features.extract$class)
t2=table(features.extract$class)

#Stratify samples according to NTP prediction: BH<0.05
prognosis <- rep(0,nrow(pred.summary)) 
for(i in seq_along(prognosis)){
  if(pred.summary[i,2] == 1 & pred.summary[i,6] < 0.05){ prognosis[i] <- "poor"}
  if(pred.summary[i,2] == 2 & pred.summary[i,6] < 0.05){ prognosis[i] <- "good"}
  if(                         pred.summary[i,6] >= 0.05){ prognosis[i] <- "intermediate"}
}

# Create Kaplan Meier curves
#read clinical data		  
meta.data <- read.table(file = input.meta.data, header=T, sep="\t", row.names=1)
km.df <- cbind(select(meta.data, c("status","time")), prognosis)

#set event and time 
set.time  = "time"    #column name of time
set.event = "status"        #column name of event

#create survival objects
poor.int.good        <- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ prognosis, data = km.df)
poor.good.remove.int <- survfit(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ prognosis, data = km.df[which(km.df$prognosis != "intermediate"),])

p.pairwise 			<- pairwise_survdiff(Surv(as.numeric(get(set.time)), as.numeric(get(set.event))) ~ prognosis, data = km.df, p.adjust.method = "none")


#plot event curve
event.plot.1 <- ggsurvplot(poor.int.good,         fun= "event", pval = surv_pvalue(poor.int.good)$pval,        conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata", linetype = "strata", ggtheme = theme_bw(),  palette = c("mediumblue", "darkgrey", "tomato")) #poor VS intermediate VS good 
event.plot.1.no.pvalue <- ggsurvplot(poor.int.good,         fun= "event",        conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata", linetype = "strata", ggtheme = theme_bw(),  palette = c("mediumblue", "darkgrey", "tomato")) #poor VS intermediate VS good 
event.plot.2 <- ggsurvplot(poor.good.remove.int,  fun= "event", pval = surv_pvalue(poor.good.remove.int)$pval, conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata", linetype = "strata", ggtheme = theme_bw(), palette = c("mediumblue", "tomato")) #poor VS good (remove intermediates)

p2=arrange_ggsurvplots(list(event.plot.1,event.plot.2), ncol=2, title= "Validation_NTP prediction")		  
p2

#create truncated Kaplan Meier curve when patients at risk drops below 10% for visualization purpose (not for p-value!)
#store p-values
pval.plot.1 <- surv_pvalue(poor.int.good)$pval
pval.plot.2	<- surv_pvalue(poor.good.remove.int)$pval

#settings
cutoff.at.risk.percentage = 10 #set percentage to cut eg. 'cutoff.at.risk = 10' for 10%
#set.event = "status" #set event to find 10% at risk

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
event.plot 				<- ggsurvplot(poor.int.good,  fun= "event", pval = surv_pvalue(poor.int.good)$pval, pval.coord = c(0,1), conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata",  palette = c("mediumblue", "darkgrey", "tomato"), xlim = c(0,cutoff.time.x.axis)) 
event.plot.no.pvalue 	<- ggsurvplot(poor.int.good,  fun= "event",        conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata",  palette = c("mediumblue", "darkgrey", "tomato"), xlim = c(0,cutoff.time.x.axis )) 
event.plot.2 			<- ggsurvplot(poor.good.remove.int,  fun= "event", pval = surv_pvalue(poor.good.remove.int)$pval, pval.coord = c(0,1), conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata", palette = c("mediumblue", "tomato"), xlim = c(0,cutoff.time.x.axis)) 
event.plot.3 			<- ggsurvplot(poor.int.good,  fun= "event", pval = paste("poor VS good \n p =", surv_pvalue(poor.good.remove.int)$pval), pval.coord = c(0,0.95), conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata", palette = c("mediumblue", "darkgrey","tomato"), xlim = c(0,cutoff.time.x.axis))

arrange_ggsurvplots(list(event.plot,event.plot.2), ncol=2, title= "TITLE HERE")	

#reset X axis with 0,5,10,15 at time axis
event.plot.3 	<- ggsurvplot(poor.int.good,  fun= "event", pval = paste("poor VS good \n p =", surv_pvalue(poor.good.remove.int)$pval), pval.coord = c(0,0.95), conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata", palette = c("mediumblue", "darkgrey","tomato"), xlim = c(0,cutoff.time.x.axis), break.time.by=5)
event.plot.3 

p=ggsurvplot(poor.int.good,
             fun= "event",
             pval = paste0("poor VS good \np=", signif(surv_pvalue(poor.good.remove.int)$pval,3)),
             pval.coord = c(0,0.80),
             pval.size=7,
             conf.int = FALSE,
             risk.table = TRUE, risk.table.fontsize=6,
             risk.table.height=0.3,
             risk.table.y.text = TRUE,
             font.x=20,
             font.y=20,
             font.tickslab=20,
             break.time.by=5,
             palette = c("mediumblue", "darkgrey","tomato"),
             xlim = c(0,cutoff.time.x.axis),
             legend = 'none',
             xlab="Time (years)",
             ylab="HCC incidence"
)
p$table <- p$table + 
  theme(axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) + ylab("Prognosis") + 
  scale_y_discrete(labels=c('poor', 'intermediate', 'good'))
p
p2 <- p
if(save.plots==TRUE){  tiff('output/plots/KM_validation.tiff', units="in", width=7, height=7, res=300, compression = 'lzw');print(p);dev.off()  }


