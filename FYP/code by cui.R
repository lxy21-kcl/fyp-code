library(dplyr)
library(tidyr)
library(lme4)
library(parameters)
library(grDevices)
library(confintr)
library(irr)
library(psych)
library(psychometric)
library(readr)
datainfo_transpose_1_20_1_ <- read_csv("C:/Users/cbh/Desktop/datainfo_transpose_1.20 (1).csv")
View(datainfo_transpose_1_20_1_)
ds<-read.csv("C:/Users/cbh/Desktop/datainfo_transpose_1.20 (1).csv",stringsAsFactors=FALSE, fileEncoding="latin1")
ds<-ds[ds$Study.Phase=="3",]

#HR 
HRds<-as.data.frame(cbind(ds$PMID,ds$Sample.size,ds$HR.for.BICR.PFS,ds$X95..CI.BICR.HR.L,ds$X95..CI.BICR.HR.U,
            ds$HR.for.INV.PFS,ds$X95..CI.INV.L,ds$X95..CI.INV.U))

#HRds<-na.omit(HRds) #if you want to use SE as weight, run this line of code

for ( i in 2:ncol(HRds)){
  HRds[,i]<-as.numeric(HRds[,i])
}
HRds$Blinded.open.label <- ds$Blinded.open.label
HRds$Cancer.Type <- ds$Cancer.Type
colnames(HRds)<-c("PMID","Sample_size","HR_BICR","HR_BICR_L","HR_BICR_U",
                  "HR_INV","HR_INV_L","HR_INV_U","Blinded.open.label", "Cancer.Type")
#"Blinded.open.label","Cancer.Type"
HRds<-HRds[!is.na(HRds$HR_BICR)&!is.na(HRds$HR_INV),] #if you want the sample size as the weight, run this line of code
#you can add other information for further analysis

#We use HR of INV/HR of BICR

#correlation
#check normality

shapiro.test(log(HRds$HR_BICR))
shapiro.test(log(HRds$HR_INV))
cor(log(HRds$HR_BICR),log(HRds$HR_INV),method="spearman") #if normality hold, use "pearson" in the method

cor.test(
  x=log(HRds$HR_BICR),
  y=log(HRds$HR_INV),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
  ) #if normality hold, use "pearson" in the method

open01_ds<-ds[ds$Blinded.open.label=="open label",]
cor.test(
  x=log(open01_ds$HR.for.BICR.PFS),
  y=log(open01_ds$HR.for.INV.PFS),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
  )#analysis open label subgroup

open02_ds<-ds[ds$Blinded.open.label=="double-blind",]
cor.test(
  x=log(open02_ds$HR.for.BICR.PFS),
  y=log(open02_ds$HR.for.INV.PFS),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
)#analysis double—blind subgroup

open03_ds<-ds[ds$Cancer.Type=="LYM",]
cor.test(
  x=log(open03_ds$HR.for.BICR.PFS),
  y=log(open03_ds$HR.for.INV.PFS),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
)#analysis Lymphoma subgroup

open04_ds<-ds[ds$Cancer.Type=="MM",]
cor.test(
  x=log(open04_ds$HR.for.BICR.PFS),
  y=log(open04_ds$HR.for.INV.PFS),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
)#analysis Myeloma subgroup

open05_ds<-ds[ds$Cancer.Type=="LEU",]
cor.test(
  x=log(open05_ds$HR.for.BICR.PFS),
  y=log(open05_ds$HR.for.INV.PFS),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
)#analysis Leukemia subgroup

open06_ds<-ds[ds$Cancer.Type=="LEU/LYM",]
cor.test(
  x=log(open06_ds$HR.for.BICR.PFS),
  y=log(open06_ds$HR.for.INV.PFS),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
)#analysis Lymphoma/Leukemia subgroup

ci_cor(log(ORds$OR_BICR),log(ORds$OR_INV),method="spearman",,type="bootstrap") #only looked at the CI

#weighted linear regression
lm<-lm(log(HR_BICR)~log(HR_INV),weights = Sample_size,data=HRds)
s1<-summary(lm)
s1$r.squared
CI.Rsq(s1$r.squared,n=nrow(HRds),k=1)

#weighted linear regression about open label subgroup
lm<-lm(log(HR.for.BICR.PFS)~log(HR.for.INV.PFS),weights = Sample.size,data=open0_ds)
s1<-summary(lm)
s1$r.squared
CI.Rsq(s1$r.squared,n=nrow(HRds),k=1)

#weighted linear regression about double-blind subgroup
lm<-lm(log(HR.for.BICR.PFS)~log(HR.for.INV.PFS),weights = Sample.size,data=open1_ds)
s1<-summary(lm)
s1$r.squared
CI.Rsq(s1$r.squared,n=nrow(HRds),k=1)

#weighted linear regression about Lymphoma subgroup
lm<-lm(log(HR.for.BICR.PFS)~log(HR.for.INV.PFS),weights = Sample.size,data=open2_ds)
s1<-summary(lm)
s1$r.squared
CI.Rsq(s1$r.squared,n=nrow(HRds),k=1)

#weighted linear regression about Myeloma subgroup
lm<-lm(log(HR.for.BICR.PFS)~log(HR.for.INV.PFS),weights = Sample.size,data=open3_ds)
s1<-summary(lm)
s1$r.squared
CI.Rsq(s1$r.squared,n=nrow(HRds),k=1)

#weighted linear regression about Leukemia subgroup
lm<-lm(log(HR.for.BICR.PFS)~log(HR.for.INV.PFS),weights = Sample.size,data=open4_ds)
s1<-summary(lm)
s1$r.squared
CI.Rsq(s1$r.squared,n=nrow(HRds),k=1)

#random effect model to calculate HR_INV/HR_BICR

#calculate SE of log HR
HRds$SE_BICR<-(log(HRds$HR_BICR_U)-log(HRds$HR_BICR_L))/(2*1.96)
HRds$SE_INV<-(log(HRds$HR_INV_U)-log(HRds$HR_INV_L))/(2*1.96)



#format the data 
BICR<-HRds[,c(1,2,3,4,5,9,10,12)]
BICR$BI<-rep("BICR",nrow(BICR))
colnames(BICR)<-c("PMID","Sample_size","HR","HR_L","HR_U","Blinded.open.label", "Cancer.Type","SE","BI")
INV<-HRds[,c(1,2,6,7,8,9,10,11)]
INV$BI<-rep("INV",nrow(INV))
colnames(INV)<-c("PMID","Sample_size","HR","HR_L","HR_U","Blinded.open.label", "Cancer.Type","SE","BI")
longHRds<-rbind(BICR,INV)

#random-effect model:
randomfit<-lmer(log(HR) ~ factor(BI)+ (1 | PMID), data = longHRds,weights=Sample_size) #OR change to SE on weights
#HR_INV/HR_BICR
s<-summary(randomfit)

exp(s$coefficients[2])
#95% CI
ci<-ci(randomfit)
exp(ci$CI_low[2])
exp(ci$CI_high[2])

#random-effect model of open label:
subdata1_ds<-longHRds[longHRds$Blinded.open.label=="open label",]
randomfit1 <- lmer(log(HR) ~ factor(BI) + (1 | PMID), data = subdata1_ds, weights = Sample_size) #OR change to SE on weights
#HR_INV/HR_BICR
s1<-summary(randomfit1)
s1$coefficients
exp(s1$coefficients[2])
#95% CI
ci<-ci(randomfit1)
exp(ci$CI_low[2])
exp(ci$CI_high[2])

#random-effect model of double-blind:
subdata2_ds<-longHRds[longHRds$Blinded.open.label=="double-blind",]
randomfit2 <- lmer(log(HR) ~ factor(BI) + (1 | PMID), data = subdata2_ds, weights = Sample_size) #OR change to SE on weights
#HR_INV/HR_BICR
s2<-summary(randomfit2)
s2$coefficients
exp(s2$coefficients[2])
#95% CI
ci<-ci(randomfit2)
exp(ci$CI_low[2])
exp(ci$CI_high[2])

#you can run the analysis for open lable/double blinded; cancer types, etc. 
#See Roche's paper for reference 



#consistency of significance
HRsig<-as.data.frame(cbind(ds$PMID,ds$Sample.size,ds$HR.INV.p.value..log.rank.test.,ds$HR.BICR.p.value))
HRsig<-na.omit(HRsig) 
for ( i in 2:ncol(HRsig)){
  HRsig[,i]<-as.numeric(HRsig[,i])
}
colnames(HRsig)<-c("PMID","Sample_size","pINV","pBICR")
HRsig$sigBICR<-ifelse(HRsig$pBICR<0.05,1,0)
HRsig$sigINV<-ifelse(HRsig$pINV<0.05,1,0)
x<-table(HRsig$sigBICR,HRsig$sigINV)
p0<-(x[1,1]+x[2,2])/nrow(HRsig)
chance1<-(sum(x[1,])/nrow(HRsig))*(sum(x[,1])/nrow(HRsig))
chance2<-(sum(x[2,])/nrow(HRsig))*(sum(x[,2])/nrow(HRsig))
pe<-chance1+chance2
kappa<-(p0-pe)/(1-pe)
kappa

rater1<-HRsig$sigBICR
rater2<-HRsig$sigINV

kappa_result <- cohen.kappa(data.frame(rater1, rater2))
print(kappa_result)

#Plot
#This part can also be used to evaluate consistency
HRds$HR_BICR_sign<-ifelse(HRds$HR_BICR_L<1&HRds$HR_BICR_U>1,0,1)
HRds$HR_INV_sign<-ifelse(HRds$HR_INV_L<1&HRds$HR_INV_U>1,0,1)
HRds<-HRds[order(HRds$HR_INV),]



x11(height=10,width=15)
op <- par(font = 1, cex = 0.9)
par(mar=c(7.5,4,4,4),mgp=c(3,1,0))
par(oma=c(2,2,2,2))

png(file="plot.png",width=1800, height=1000,res=144)
plot(c(0,nrow(HRds)), rep(1,2), log="y", ylim=c(0.05,3), xlim=c(0,nrow(HRds)+1), axes=FALSE, xlab="Trials", ylab="HR estimate (95%CI)", type="l",cex.lab=1.2)
axis(1, at = 0:(nrow(HRds)-1),labels = c(1:(nrow(HRds))), las = 2, pos=c(0.05,-1))
axis(2, at = c(0.05,0.1,0.2,0.5,1,2,3), labels = c( 0.05,0.1,0.2,0.5,1,2,3), pos=c(-0.7,0))
legend(x=4,y=3, legend=c("INV", "BICR"), pch=c(1,19),cex = 0.8)

for(i in 1:nrow(HRds)){
  if (HRds$HR_INV_sign[i] == 0){
    points(i-1-0.2, (HRds$HR_INV[i]),cex=0.7,lwd=1)
    segments(i-1-0.2, (HRds$HR_INV_L[i]), i-1-0.2, (HRds$HR_INV_U[i]), col = gray(0.5))
  }
  if (HRds$HR_INV_sign[i] == 1){
    points(i-1-0.2, (HRds$HR_INV[i]),cex=0.7,lwd=1)
    segments(i-1-0.2, (HRds$HR_INV_L[i]), i-1-0.2, (HRds$HR_INV_U[i]), col="#0072B2")
  }
}

for(i in 1:nrow(HRds)){
  if (HRds$HR_BICR_sign[i] == 0){
    points(i-1+0.2, (HRds$HR_BICR[i]), pch=16,cex=0.7)
    segments(i-1+0.2, (HRds$HR_BICR_L[i]), i-1+0.2, (HRds$HR_BICR_U[i]), col = gray(0.5))
  }
  if (HRds$HR_BICR_sign[i] == 1){
    points(i-1+0.2, (HRds$HR_BICR[i]), pch=16,cex=0.7)
    segments(i-1+0.2, (HRds$HR_BICR_L[i]), i-1+0.2, (HRds$HR_BICR_U[i]), col = "#E69F00")
  }
}
dev.off()



#OR
datainfo_transpose_1_20_1_ <- read_csv("C:/Users/cbh/Desktop/datainfo_transpose_1.20 (1).csv")
View(datainfo_transpose_1_20_1_)
ds<-read.csv("C:/Users/cbh/Desktop/datainfo_transpose_1.20 (1).csv",stringsAsFactors=FALSE, fileEncoding="latin1")
ds<-ds[ds$Study.Phase=="3",]

ORds<-as.data.frame(cbind(ds$PMID,ds$Sample.size,ds$BICR.ORR.experimental,ds$X95..CI.BICR.ORR.L.experimental,
                          ds$X95..CI.BICR.ORR.U.experimental,ds$BICR.ORR.control,ds$X95..CI.BICR.ORR.L.control,ds$X95..CI.BICR.ORR.U.control,
                          ds$INV.ORR.experimental,ds$X95..CI.INV.ORR.L.experimental,ds$X95..CI.INV.ORR.U.experimental,
                          ds$INV.ORR.control,ds$X95..CI.INV.ORR.L.control,
                          ds$X95..CI.INV.ORR.U.control))

for ( i in 2:ncol(ORds)){
  ORds[,i]<-as.numeric(ORds[,i])
}
ORds$Blinded.open.label<-ds$Blinded.open.label
ORds$Cancer.Type<-ds$Cancer.Type
colnames(ORds)<-c("PMID","Sample_size","OR_BICR_exp","OR_BICR_L_exp","OR_BICR_U_exp",
                  "OR_BICR_c","OR_BICR_L_c","OR_BICR_U_c",
                  "OR_INV_exp","OR_INV_L_exp","OR_INV_U_exp",
                  "OR_INV_c","OR_INV_L_c","OR_INV_U_c","Blinded.open.label","Cancer.Type")

ORds$OR_BICR<-(ORds$OR_BICR_exp*(1-ORds$OR_BICR_exp))/(ORds$OR_BICR_c*(1-ORds$OR_BICR_c))
ORds$OR_INV<-(ORds$OR_INV_exp*(1-ORds$OR_INV_exp))/(ORds$OR_INV_c*(1-ORds$OR_INV_c))
ORds<-ORds[!is.na(ORds$OR_BICR)&!is.na(ORds$OR_INV),]
#you can add other information for further analysis

#We use OR of INV/OR of BICR

#correlation
#check normality

shapiro.test(log(ORds$OR_BICR))
shapiro.test(log(ORds$OR_INV))
cor(log(ORds$OR_BICR),log(ORds$OR_INV))
ci_cor(log(ORds$OR_BICR),log(ORds$OR_INV),method="pearson")

open11_ds<-ORds[ORds$Blinded.open.label=="open label",]
cor.test(
  x=log(open11_ds$OR_BICR),
  y=log(open11_ds$OR_INV),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
)#analysis open label subgroup

open12_ds<-ORds[ORds$Blinded.open.label=="double-blind",]
cor.test(
  x=log(open12_ds$OR_BICR),
  y=log(open12_ds$OR_INV),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
)#analysis double-blind subgroup

open13_ds<-ORds[ORds$Cancer.Type=="LYM",]
cor.test(
  x=log(open13_ds$OR_BICR),
  y=log(open13_ds$OR_INV),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
)#analysis Lymphoma subgroup

open14_ds<-ORds[ORds$Cancer.Type=="MM",]
cor.test(
  x=log(open14_ds$OR_BICR),
  y=log(open14_ds$OR_INV),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
)#analysis Myeloma subgroup

open15_ds<-ORds[ORds$Cancer.Type=="LEU",]
cor.test(
  x=log(open15_ds$OR_BICR),
  y=log(open15_ds$OR_INV),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
)#analysis Leukemia subgroup

open16_ds<-ORds[ORds$Cancer.Type=="LEU/LYM",]
cor.test(
  x=log(open16_ds$OR_BICR),
  y=log(open16_ds$OR_INV),
  method="pearson",
  exact=FALSE,
  conf.level=0.95
)#analysis Lymphoma/Leukemia subgroup

#weighted linear regression
lm<-lm(log(OR_BICR)~log(OR_INV),weights = Sample_size,data=ORds)
s1<-summary(lm)
s1$r.squared
CI.Rsq(s1$r.squared,n=nrow(HRds),k=1)

#weighted linear regression about open label subgroup
lm<-lm(log(OR_BICR)~log(OR_INV),weights = Sample_size,data=open11_ds)
s1<-summary(lm)
s1$r.squared
CI.Rsq(s1$r.squared,n=nrow(HRds),k=1)

#weighted linear regression about double-blind subgroup
lm<-lm(log(OR_BICR)~log(OR_INV),weights = Sample_size,data=open12_ds)
s1<-summary(lm)
s1$r.squared
CI.Rsq(s1$r.squared,n=nrow(HRds),k=1)

#weighted linear regression about Lymphoma subgroup
lm<-lm(log(OR_BICR)~log(OR_INV),weights = Sample_size,data=open13_ds)
s1<-summary(lm)
s1$r.squared
CI.Rsq(s1$r.squared,n=nrow(HRds),k=1)

#weighted linear regression about Myeloma subgroup
lm<-lm(log(OR_BICR)~log(OR_INV),weights = Sample_size,data=open14_ds)
s1<-summary(lm)
s1$r.squared
CI.Rsq(s1$r.squared,n=nrow(HRds),k=1)

#weighted linear regression about Leukemia subgroup
lm<-lm(log(OR_BICR)~log(OR_INV),weights = Sample_size,data=open15_ds)
s1<-summary(lm)
s1$r.squared
CI.Rsq(s1$r.squared,n=nrow(HRds),k=1)


#random effect model to calculate OR_Bicr/OR_LE


#format the data 
BICR<-ORds[,c(1,2,15)]
BICR$BI<-rep("BICR",nrow(BICR))
colnames(BICR)<-c("PMID","Sample_size","OR","BI")
INV<-ORds[,c(1,2,16)]
INV$BI<-rep("INV",nrow(INV))
colnames(INV)<-c("PMID","Sample_size","OR","BI")
longORds<-rbind(BICR,INV)

#random-effect model:
randomfit2<-lmer(log(OR) ~ factor(BI)+ (1 | PMID), data = longORds,weights=Sample_size)
#OR_INV/OR_BICR
s2<-summary(randomfit2)
s2
exp(s2$coefficients[2])
#95% CI
ci<-ci(randomfit2)
exp(ci$CI_low[2])
exp(ci$CI_high[2])

#you can run the analysis for open lable/double blinded; cancer types, etc. 
#See Roche's paper for reference 



