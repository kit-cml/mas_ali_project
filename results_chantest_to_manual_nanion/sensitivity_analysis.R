rm(list=ls())
library(rms)
set.seed(1)
library(ggplot2)
source("functions.R")
#=================
#Input information

#--- specify command line arguments
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, c("-d", "--datadir"), default = "data/metrics_manual_tms.csv", help = "Filepath to metrics file")
parser <- add_option(parser, c("-m", "--dataset_name"), default = "manual", help = "Name of dataset")
parser <- add_option(parser, c("-o", "--outdir"), default = "results_sensitivity", help = "Output directory")
parser <- add_option(parser, c("-i", "--influential_drugs"), default = 10, type="integer", help = "Number of influential drugs to be captured for each TMS threshold analysis")
parser <- add_option(parser, c("-c", "--cpr"), default = 0.75, type="double", help = "Correct prediciton rate")
args <- parse_args(parser)

datadir <- args$datadir
dataset_name <- args$dataset_name
outdir <- args$outdir
influential_drugs <- args$influential_drugs
cpr <- args$cpr

sprintf("datadir = %s",datadir)
sprintf("dataset_name = %s", dataset_name)
sprintf("outdir = %s",outdir)
sprintf("influential_drugs = %s",as.character(influential_drugs))
sprintf("cpr = %s",as.character(cpr))

# Data directory location
#datadir <- "data/metrics_nanion_scaled_w_tms_dvdtmax-APD90-camax-carest.csv"
#dataset_name <- "nanion"

# Perturbation percentage 
variations <- seq(0.75,1.25,by=0.05)

#outdir<-sprintf("results_sensitivity")
#influential_drugs <- 10
#cpr <- 0.75

#=================

scale<-"sara"
CL<-2000
metric<-"tms"
dataset<-"chantest"
last_dose<-4
dose<-last_dose
fitmethod<-"lrm"
#fitmethod<-"ordLORgee"
#fitmethod<-"repolr"

colorvec<-c("red","dark green","blue","black","cyan","purple","brown")
colordf<-data.frame(
  drug = c("dofetilide","bepridil","quinidine","sotalol","ibutilide", "azimilide", "disopyramide", "vandetanib",
           "cisapride","terfenadine","ondansetron","chlorpromazine","clarithromycin","risperidone","domperidone","astemizole", "pimozide","droperidol","clozapine",
           "ranolazine","mexiletine","diltiazem","verapamil","metoprolol","loratadine","tamoxifen","nifedipine","nitrendipine"),
  classidx=c(rep(2,8),rep(1,11),rep(0,9)),
  coloridx=c(rep(1,8),rep(3,11),rep(2,9)),
  isTraining=c(1,1,1,1,0,0,0,0,
               1,1,1,1,0,0,0,0,0,0,0,
               1,1,1,1,0,0,0,0,0)
)

system(paste0("mkdir -p ",outdir))
figdir<-sub("results","figs",outdir)
system(paste0("mkdir -p ",figdir))

if(fitmethod=="polr"){
  library(MASS)
  try_fit<-try_polr
}else if(fitmethod=="lrm"){
  library(rms)
  try_fit<-try_lrm
}else if(fitmethod=="repolr"){
  require(repolr)
  try_fit<-try_repolr
}else if(fitmethod=="ordLORgee"){
  require(multgee)
  try_fit<-try_ordLORgee
}

# read in dataset

thrdf<-data.frame()

if(metric !="model5"){
  df<- as.data.frame(read.csv(datadir))  
  #mei head
  #df<-df[df$drug!="disopyramide"&df$drug!="domperidone",]
  #mei over
  df<-df[df$drug!="control" & df$dose>0,]
  todrop<-is.na(df$max_dv)#df$noDepols==250
  print(sprintf("Number of depolarization failures: %d/%d",sum(todrop),nrow(df)))
  df<-df[!todrop,]
  df$class<-df$max_dv>0
  # df$qNet <- df$qNet/1000
  print(sprintf("Number of EADs: %d/%d",sum(df$class),nrow(df)))
  
  if(dataset=="original"){
    nboots<-0
  }else{
    #nboots<-max(df$boot) # using all bootstraps
    nboots<-2000
  }
  
  cols<-c("drug","sample")
  
  # create thresholds for each range of doses
  print(sprintf("CL = %d ms, nboot = %d, dose <= %d Cmax",CL,nboots,dose))
  # print(sprintf("%s, %s scale, CL = %d ms, nboot = %d, dose <= %d Cmax",dataset,scale,CL,nboots,dose))
  
  brows<-df$dose<=dose
  wtdf<-aggregate(df[brows,metric,drop=FALSE], by=as.list(df[brows,cols]), mean) # equally weighted across doses
  wtdf$drug<-as.factor(wtdf$drug)
}else{                            #if metric == model5
  if(dataset=="pureHTS"){                   #shouldn't use "hybrid" as dataset here!
    Hilldataset<- "pureHTS"
  }else if(dataset=="pureManual"){
    Hilldataset<-"pureManual"          #note I didn't use validationManual_Hill_fitting here
  }
  wtdf<-data.frame(drug=rep(colordf$drug, each=2000),sample=rep(1:2000,28),model5=NA)
  for(drug in colordf$drug){
    IC50table<- read.delim(paste0("/scratch/lizhi/cardio/UQ/",Hilldataset,"_Hill_fitting/results/",drug,"/IC50_samples.csv"),sep=",")
    wtdf$model5[wtdf$drug==drug]<- with(IC50table, -log(ICaL_IC50) - -log(hERG_IC50))
  }#for drug
}#if metric !=model5
# get 95% CI
q2.5<-aggregate(wtdf[,metric,drop=FALSE], by=list(drug=wtdf[,c("drug")]),
                FUN=function(x) quantile(x,probs=0.025,na.rm=TRUE))
q97.5<-aggregate(wtdf[,metric,drop=FALSE], by=list(drug=wtdf[,c("drug")]),
                 FUN=function(x) quantile(x,probs=0.975,na.rm=TRUE))

q50<-aggregate(wtdf[,metric,drop=FALSE], by=list(drug=wtdf[,c("drug")]),
               FUN=function(x) quantile(x,probs=0.5,na.rm=TRUE))                                                  
cidf<-merge(q2.5,q97.5,by=c("drug"),suffixes=c("_0.025","_0.975"))
cidf<-merge(cidf, q50, by="drug")
colnames(cidf)[4]<-paste0(metric,"_0.5")
outdf<-cidf

outfile<-sprintf("%s/%s_allCIthresholds.rds",outdir,metric)
saveRDS(outdf, outfile)


outdf<- merge(colordf, outdf,by="drug")
cidf<-outdf

desiredorder<-drugnames<-c("ibutilide","vandetanib","bepridil","azimilide","dofetilide","quinidine","sotalol",  "disopyramide", 
                           "domperidone","pimozide","droperidol","cisapride","terfenadine","astemizole","ondansetron","clarithromycin","chlorpromazine","clozapine","risperidone", 
                           "tamoxifen","loratadine","verapamil","metoprolol","nitrendipine","diltiazem","ranolazine","nifedipine","mexiletine")

outdf$drug<- factor(outdf$drug, levels=rev(desiredorder))
outdf$class<- 2-outdf$classidx
lower<-paste(metric,0.025,sep="_")
upper<-paste(metric,0.975,sep="_")
middle<-paste(metric,0.5,sep="_")
lower_bds<-aggregate(outdf[,lower], by=list(risk=outdf$class), FUN=min)
upper_bds<-aggregate(outdf[,upper], by=list(risk=outdf$class), FUN=max)
lower_bds<-lower_bds[match(c(0:2),lower_bds$risk),]
upper_bds<-upper_bds[match(c(0:2),upper_bds$risk),]
thresholds<-apply(cbind(head(upper_bds$x,-1),tail(lower_bds$x,-1)), 1, max)
thresholds<-sort(thresholds)
names(thresholds)<-c("High","Intermediate")

newrow<-data.frame(last_dose=dose, type="95CI", t(thresholds))
thrdf<-rbind(thrdf,newrow)

figfile<-sub("results","figs",sub(".rds",".pdf",outfile))
pdf(figfile, width=8, height=4)
p<-ggplot(outdf, aes(x=drug, color=as.character(class)))    #ggplot2 doesn't like color idx being numbers
p<-p+geom_errorbar(aes_string(ymin=lower, ymax=upper))
for(thresh in thresholds)
  p<-p+geom_hline(yintercept=thresh, linetype="dotted")
p<-p+ylab(sprintf("%s_score_1-%gX_Cmax",metric,dose))
p<-p+coord_flip()
p<-p+scale_color_brewer(NULL, palette="Set1")
p<-p+theme_bw()
print(p)
dev.off()

# get ordinal logistic regression thresholds
wtdf<-merge(wtdf,colordf)

wtdf$class <- ordered(2-wtdf$classidx)  #change the order because in Kelly's code "high" is first class (0).
#but in mycompute_TdP_error.R "high" is 2 as planned

# fit all data using single predictor
datadf<-wtdf[,c("drug","class","sample",metric)]
#datadf$drug<- factor(datadf$drug, levels=rev(desiredorder))
maxsamp <- 2000
datadf<-datadf[datadf$sample<=maxsamp,]
maxreplacedidx<-0        #to break the perfect correlation within drugs for GEE methods
#for ordinary lrm can use 0
if(maxreplacedidx !=0){                                               
  test<-do.call(rbind,by(datadf, datadf$drug, function(x) {
    class<-x$class; fullvec<-0:2;
    for(i in 1:maxreplacedidx){
      c<-class[i]; idx<-fullvec%in%c;class[i]<-sample(fullvec[!idx],1);
    }
    x$class<-class;return(x)
  }))
  datadf<-test[order(test$drug,test$sample),]  #sort by sample (time) is required by repolr but not ordLORgee
}


if(fitmethod=="lrm"){
  dd<-datadist(datadf)
  options(datadist="dd")
}else if(fitmethod=="repolr"){             #repolr doesn't like factors! and doesn't like class being 0?
  
  datadf$class<-as.integer(datadf$class)  #if -1 then the values are the same as factor (0,1,2)
}
lmod<-try_fit(datadf, metric, 0) # no penalty
print(lmod)

if(inherits(lmod, "try-error")){
  print(sprintf("fitting dose %d failed! skipping...",dose))
  next
}

# save coefficients
if(fitmethod=="polr"){
  print(sprintf("Convergence code: %d, Number of iterations: %d",lmod$convergence,lmod$niter))
  cf<-lmod$coefficients[[metric]]
  ints<-lmod$zeta
  cfvec<-c()
  for(kint in 1:length(ints))
    cfvec[[paste0("intercept",kint)]]<-ints[[kint]]
  cfvec[["slope"]]<-cf
}else{
  print(sprintf("Convergence failure: %s",lmod$fail))
  cf<-coefficients(lmod)
  cfvec<-c()
  for(kint in 1:(length(cf)-1))
    cfvec[[paste0("intercept",kint)]]<-cf[[kint]]
  cfvec[["slope"]]<-cf[[length(cf)]]
}
t1<- -cfvec[["intercept1"]]/cfvec[["slope"]]
t2<- -cfvec[["intercept2"]]/cfvec[["slope"]]
#use formal math
ebeta1<-exp(-cfvec[["intercept1"]]); ebeta2<- exp(-cfvec[["intercept2"]])
t1 <- log(ebeta1*ebeta2/-(2*ebeta1-ebeta2))/cfvec[["slope"]]
t2 <- log(exp(-cfvec[["intercept2"]])-2*exp(-cfvec[["intercept1"]]))/cfvec[["slope"]]

thresholds<-c(High=t2, Intermediate=t1)
print(thresholds)

high1<-t2;intermediate1<-t1;

# Check the wrongly predicted drugs by looking at the median values
# of prediction results

y0 <- datadf$class

probs <- matrix(predict(lmod, newdata = datadf, type = "fitted.ind"), ncol = length(levels(y0)))

#sometimes APD90 is NA, giving rise to NA in probs
idx <- apply(probs, 1, function(x) any(is.na(x)))
probs[idx, ] <- matrix(rep(c(0, 0, 1), sum(idx)), nrow = sum(idx), byrow = T) #this is to assume all NAs are due to repolarization failure, and thus high risk

yPred <- apply(probs, 1, function(x) which.max(x)) - 1
pred_err <- yPred == y0

drug_predict <- data.frame(
  drug = as.character(datadf$drug),
  pred = yPred,
  label = datadf$class,
  pred_err = pred_err
)

selected_drugs <- data.frame()
for (drug in unique(drug_predict$drug)) {
  temp <- drug_predict$pred_err[drug_predict$drug==drug]
  correct_pred_rate <- sum(temp)/length(temp)
  if (correct_pred_rate >= cpr) {
    temp_df <- data.frame(
      drug = drug,
      cpr = correct_pred_rate
    )
    selected_drugs <- rbind(selected_drugs,temp_df)
  }
}

xmdf<- as.data.frame(read.csv(datadir))

xmthreshold1<-list()#28drug
xmthreshold2<-list()
xmdrug<-unique(xmdf$drug)

for (xmn in 1:length(xmdrug)){
  thresholdslist1<-list()#16thresholds
  thresholdslist2<-list()
  xmnum<-variations 
  
  for (m in 1:length(xmnum)){
    
    xm1<-xmdf[xmdf$drug!=xmdrug[xmn],]
    xm2<-xmdf[xmdf$drug==xmdrug[xmn],]
    xm2[paste0(metric)]<-xm2[paste0(metric)]*xmnum[m]
    # xm2$tms<-xm2$tms*xmnum[m]
    # xm2$qNet<-xm2$qNet*xmnum[m] # what is this for?
    #meio
    xmwtdf<-rbind(xm1,xm2)
    
    
    if(fitmethod=="polr"){
      library(MASS)
      try_fit<-try_polr
    }else if(fitmethod=="lrm"){
      library(rms)
      try_fit<-try_lrm
    }else if(fitmethod=="repolr"){
      require(repolr)
      try_fit<-try_repolr
    }else if(fitmethod=="ordLORgee"){
      require(multgee)
      try_fit<-try_ordLORgee
    }
    
    # read in dataset
    
    thrdf<-data.frame()
    
    if(metric !="model5"){
      df<-xmwtdf  
      df<-df[df$drug!="control" & df$dose>0,]
      todrop<-is.na(df$max_dv)#df$noDepols==250
      print(sprintf("Number of depolarization failures: %d/%d",sum(todrop),nrow(df)))
      df<-df[!todrop,]
      df$class<-df$max_dv>0 ### SAMPAI DI SINI!!!
      # df$qNet <- df$qNet/1000
      print(sprintf("Number of EADs: %d/%d",sum(df$class),nrow(df)))
      
      
      if(dataset=="original"){
        nboots<-0
      }else{
        #nboots<-max(df$boot) # using all bootstraps
        nboots<-2000
      }
      
      cols<-c("drug","sample")
      
      
      # create thresholds for each range of doses
      print(sprintf("%s, %s scale, CL = %d ms, nboot = %d, dose <= %d Cmax",dataset,scale,CL,nboots,dose))
      
      brows<-df$dose<=dose
      wtdf<-aggregate(df[brows,metric,drop=FALSE], by=as.list(df[brows,cols]), mean) # equally weighted across doses
      wtdf$drug<-as.factor(wtdf$drug)
    }else{                            #if metric == model5
      if(dataset=="pureHTS"){                   #shouldn't use "hybrid" as dataset here!
        Hilldataset<- "pureHTS"
      }else if(dataset=="pureManual"){
        Hilldataset<-"pureManual"          #note I didn't use validationManual_Hill_fitting here
      }
      wtdf<-data.frame(drug=rep(colordf$drug, each=2000),sample=rep(1:2000,28),model5=NA)
      for(drug in colordf$drug){
        IC50table<- read.delim(paste0("/scratch/lizhi/cardio/UQ/",Hilldataset,"_Hill_fitting/results/",drug,"/IC50_samples.csv"),sep=",")
        wtdf$model5[wtdf$drug==drug]<- with(IC50table, -log(ICaL_IC50) - -log(hERG_IC50))
      }#for drug
      
      
      
      
    }#if metric !=model5
    # get 95% CI
    q2.5<-aggregate(wtdf[,metric,drop=FALSE], by=list(drug=wtdf[,c("drug")]),
                    FUN=function(x) quantile(x,probs=0.025,na.rm=TRUE))
    q97.5<-aggregate(wtdf[,metric,drop=FALSE], by=list(drug=wtdf[,c("drug")]),
                     FUN=function(x) quantile(x,probs=0.975,na.rm=TRUE))
    
    q50<-aggregate(wtdf[,metric,drop=FALSE], by=list(drug=wtdf[,c("drug")]),
                   FUN=function(x) quantile(x,probs=0.5,na.rm=TRUE))                                                  
    cidf<-merge(q2.5,q97.5,by=c("drug"),suffixes=c("_0.025","_0.975"))
    cidf<-merge(cidf, q50, by="drug")
    colnames(cidf)[4]<-paste0(metric,"_0.5")
    outdf<-cidf
    
    outfile<-sprintf("%s/%s_allCIthresholds.rds",outdir,metric)
    saveRDS(outdf, outfile)
    
    
    outdf<- merge(colordf, outdf,by="drug")
    cidf<-outdf
    
    desiredorder<-drugnames<-c("ibutilide","vandetanib","bepridil","azimilide","dofetilide","quinidine","sotalol",  "disopyramide", 
                               "domperidone","pimozide","droperidol","cisapride","terfenadine","astemizole","ondansetron","clarithromycin","chlorpromazine","clozapine","risperidone", 
                               "tamoxifen","loratadine","verapamil","metoprolol","nitrendipine","diltiazem","ranolazine","nifedipine","mexiletine")
    
    
    outdf$drug<- factor(outdf$drug, levels=rev(desiredorder))
    outdf$class<- 2-outdf$classidx
    lower<-paste(metric,0.025,sep="_")
    upper<-paste(metric,0.975,sep="_")
    middle<-paste(metric,0.5,sep="_")
    lower_bds<-aggregate(outdf[,lower], by=list(risk=outdf$class), FUN=min)
    upper_bds<-aggregate(outdf[,upper], by=list(risk=outdf$class), FUN=max)
    lower_bds<-lower_bds[match(c(0:2),lower_bds$risk),]
    upper_bds<-upper_bds[match(c(0:2),upper_bds$risk),]
    thresholds<-apply(cbind(head(upper_bds$x,-1),tail(lower_bds$x,-1)), 1, max)
    thresholds<-sort(thresholds)
    names(thresholds)<-c("High","Intermediate")
    
    newrow<-data.frame(last_dose=dose, type="95CI", t(thresholds))
    thrdf<-rbind(thrdf,newrow)
    
    figfile<-sub("results","figs",sub(".rds",".pdf",outfile))
    pdf(figfile, width=8, height=4)
    p<-ggplot(outdf, aes(x=drug, color=as.character(class)))    #ggplot2 doesn't like color idx being numbers
    p<-p+geom_errorbar(aes_string(ymin=lower, ymax=upper))
    for(thresh in thresholds)
      p<-p+geom_hline(yintercept=thresh, linetype="dotted")
    p<-p+ylab(sprintf("%s_score_1-%gX_Cmax",metric,dose))
    p<-p+coord_flip()
    p<-p+scale_color_brewer(NULL, palette="Set1")
    p<-p+theme_bw()
    print(p)
    dev.off()
    
    # get ordinal logistic regression thresholds
    wtdf<-merge(wtdf,colordf)
    
    wtdf$class <- ordered(2-wtdf$classidx)  #change the order because in Kelly's code "high" is first class (0).
    #but in mycompute_TdP_error.R "high" is 2 as planned
    
    # fit all data using single predictor
    datadf<-wtdf[,c("drug","class","sample",metric)]
    #datadf$drug<- factor(datadf$drug, levels=rev(desiredorder))
    maxsamp <- 2000
    datadf<-datadf[datadf$sample<=maxsamp,]
    maxreplacedidx<-0                           #to break the perfect correlation within drugs for GEE methods
    #for ordinary lrm can use 0
    if(maxreplacedidx !=0){                                               
      test<-do.call(rbind,by(datadf, datadf$drug, function(x) {
        class<-x$class; fullvec<-0:2;
        for(i in 1:maxreplacedidx){
          c<-class[i]; idx<-fullvec%in%c;class[i]<-sample(fullvec[!idx],1);
        }
        x$class<-class;return(x)
      }))
      datadf<-test[order(test$drug,test$sample),]  #sort by sample (time) is required by repolr but not ordLORgee
    }
    
    
    if(fitmethod=="lrm"){
      dd<-datadist(datadf)
      options(datadist="dd")
    }else if(fitmethod=="repolr"){             #repolr doesn't like factors! and doesn't like class being 0?
      
      datadf$class<-as.integer(datadf$class)  #if -1 then the values are the same as factor (0,1,2)
    }
    lmod<-try_fit(datadf, metric, 0) # no penalty
    print(lmod)
    
    if(inherits(lmod, "try-error")){
      print(sprintf("fitting dose %d failed! skipping...",dose))
      next
    }
    
    # save coefficients
    if(fitmethod=="polr"){
      print(sprintf("Convergence code: %d, Number of iterations: %d",lmod$convergence,lmod$niter))
      cf<-lmod$coefficients[[metric]]
      ints<-lmod$zeta
      cfvec<-c()
      for(kint in 1:length(ints))
        cfvec[[paste0("intercept",kint)]]<-ints[[kint]]
      cfvec[["slope"]]<-cf
    }else{
      print(sprintf("Convergence failure: %s",lmod$fail))
      cf<-coefficients(lmod)
      cfvec<-c()
      for(kint in 1:(length(cf)-1))
        cfvec[[paste0("intercept",kint)]]<-cf[[kint]]
      cfvec[["slope"]]<-cf[[length(cf)]]
    }
    t1<- -cfvec[["intercept1"]]/cfvec[["slope"]]
    t2<- -cfvec[["intercept2"]]/cfvec[["slope"]]
    #use formal math
    ebeta1<-exp(-cfvec[["intercept1"]]); ebeta2<- exp(-cfvec[["intercept2"]])
    t1 <- log(ebeta1*ebeta2/-(2*ebeta1-ebeta2))/cfvec[["slope"]]
    t2 <- log(exp(-cfvec[["intercept2"]])-2*exp(-cfvec[["intercept1"]]))/cfvec[["slope"]]
    
    #thresholds<-c(High=t2, Intermediate=t1)
    
    #meih
    
    
    
    thresholdslist1[[m]]<-c(xmdrug[xmn],xmnum[m],t1,intermediate1,"threshold2")
    
    thresholdslist2[[m]]<-c(xmdrug[xmn],xmnum[m],t2,high1,"threshold1")
    print(thresholdslist1[[m]])
    
  }
  xmthreshold1[[xmn]]<-do.call(rbind,thresholdslist1)
  xmthreshold2[[xmn]]<-do.call(rbind,thresholdslist2)
}

xmthreshold11<-do.call(rbind,xmthreshold1)
xmthreshold12<-do.call(rbind,xmthreshold2)
xmthreshold13<-rbind(xmthreshold11,xmthreshold12)
write.csv(xmthreshold13,paste0(outdir,"/",dataset_name,"_28_sens.csv"))

# Convert character type to numbers
sens<-data.frame(xmthreshold13)
sens$X2<- as.numeric(sens$X2)
sens$X3 <- as.numeric(sens$X3)
sens$X4 <- as.numeric(sens$X4)
sens$X6 <- abs((sens$X3-sens$X4)/sens$X4) # Threshold changes

# Rank all drugs based on threshold changes
drug1_ranks <- data.frame(drug=subset(sens,sens$X2==0.75 & sens$X5=="threshold1")[,"X1"])
tempn <- subset(sens,sens$X2==0.75 & sens$X5=="threshold1")
tempp <- subset(sens,sens$X2==1.25 & sens$X5=="threshold1")
drug1_ranks$rank_25 <- rank(-(tempn$X6+tempp$X6)/2.0)
tempn <- subset(sens,sens$X2==0.8 & sens$X5=="threshold1")
tempp <- subset(sens,sens$X2==1.2 & sens$X5=="threshold1")
drug1_ranks$rank_20 <- rank(-(tempn$X6+tempp$X6)/2.0)
tempn <- subset(sens,sens$X2==0.85 & sens$X5=="threshold1")
tempp <- subset(sens,sens$X2==1.15 & sens$X5=="threshold1")
drug1_ranks$rank_15 <- rank(-(tempn$X6+tempp$X6)/2.0)
tempn <- subset(sens,sens$X2==0.9 & sens$X5=="threshold1")
tempp <- subset(sens,sens$X2==1.1 & sens$X5=="threshold1")
drug1_ranks$rank_10 <- rank(-(tempn$X6+tempp$X6)/2.0)
tempn <- subset(sens,sens$X2==0.95 & sens$X5=="threshold1")
tempp <- subset(sens,sens$X2==1.05 & sens$X5=="threshold1")
drug1_ranks$rank_5 <- rank(-(tempn$X6+tempp$X6)/2.0)
drug1_ranks$rank_mean <- 1/5.0 * (drug1_ranks$rank_25 + drug1_ranks$rank_20 + drug1_ranks$rank_15 + drug1_ranks$rank_10 + drug1_ranks$rank_5)
drug1_ranks$rank_final <- rank(drug1_ranks$rank_mean, ties.method = "first")

drug2_ranks <- data.frame(drug=subset(sens,sens$X2==0.75 & sens$X5=="threshold2")[,"X1"])
tempn <- subset(sens,sens$X2==0.75 & sens$X5=="threshold2")
tempp <- subset(sens,sens$X2==1.25 & sens$X5=="threshold2")
drug2_ranks$rank_25 <- rank(-(tempn$X6+tempp$X6)/2.0)
tempn <- subset(sens,sens$X2==0.8 & sens$X5=="threshold2")
tempp <- subset(sens,sens$X2==1.2 & sens$X5=="threshold2")
drug2_ranks$rank_20 <- rank(-(tempn$X6+tempp$X6)/2.0)
tempn <- subset(sens,sens$X2==0.85 & sens$X5=="threshold2")
tempp <- subset(sens,sens$X2==1.15 & sens$X5=="threshold2")
drug2_ranks$rank_15 <- rank(-(tempn$X6+tempp$X6)/2.0)
tempn <- subset(sens,sens$X2==0.9 & sens$X5=="threshold2")
tempp <- subset(sens,sens$X2==1.1 & sens$X5=="threshold2")
drug2_ranks$rank_10 <- rank(-(tempn$X6+tempp$X6)/2.0)
tempn <- subset(sens,sens$X2==0.95 & sens$X5=="threshold2")
tempp <- subset(sens,sens$X2==1.05 & sens$X5=="threshold2")
drug2_ranks$rank_5 <- rank(-(tempn$X6+tempp$X6)/2.0)
drug2_ranks$rank_mean <- 1/5.0 * (drug2_ranks$rank_25 + drug2_ranks$rank_20 + drug2_ranks$rank_15 + drug2_ranks$rank_10 + drug2_ranks$rank_5)
drug2_ranks$rank_final <- rank(drug2_ranks$rank_mean, ties.method = "first")

write.csv(drug1_ranks,paste0(outdir,"/",dataset_name,"_threshold1_ranks.csv"), row.names = FALSE)
write.csv(drug2_ranks,paste0(outdir,"/",dataset_name,"_threshold2_ranks.csv"),row.names = FALSE)

# Remove duplicates among the top drugs based on "influential-drugs"
drugs1 <- subset(drug1_ranks,drug1_ranks$rank_final<=influential_drugs)
drugs2 <- subset(drug2_ranks,drug2_ranks$rank_final<=influential_drugs)
drugs <- union(drugs1$drug,drugs2$drug)

# Remove incorrectly predicted drugs
final_drugs <- data.frame(
                drug = intersect(drugs,as.character(selected_drugs$drug))
                 )
for (row in 1:nrow(final_drugs)) {
  drug <- final_drugs$drug[row]
  final_drugs$risk[row] <- colordf$classidx[colordf$drug == drug]
}
final_drugs <- final_drugs[order(final_drugs$risk, decreasing = TRUE),]
print(final_drugs)
write.csv(final_drugs,paste0(outdir,"/",dataset_name,"_drug_candidates.csv"), row.names = FALSE)