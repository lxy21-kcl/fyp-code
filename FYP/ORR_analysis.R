library(readxl)
ds<-read_excel("datainfo_transpose_1.20.xlsx")
ORds<-as.data.frame(cbind(ds$PMID,ds$`Sample size`,ds$`BICR-ORR-experimental`,
                          ds$`95% CI BICR ORR L-experimental`,ds$`95% CI BICR ORR U-experimental`,
                          ds$`INV-ORR-experimental`,ds$`95% CI INV ORR L-experimental`,
                          ds$`95% CI INV ORR U-experimental`))

ORds$V2<-as.numeric(ORds$V2)
for ( i in 3:ncol(ORds)){
  ORds[,i]<-as.numeric(ORds[,i])/100
}

colnames(ORds)<-c("PMID","Sample_size","OR_BICR_exp","OR_BICR_L_exp","OR_BICR_U_exp",
                  "OR_INV_exp","OR_INV_L_exp","OR_INV_U_exp")

ORds<-na.omit(ORds)
ORds<-ORds[ORds$Sample_size>10,]
ORds<-ORds[ORds$Sample_size*ORds$OR_BICR_exp>5,]
# Load the metafor package
library(metafor)
install.packages("meta")   
library(meta) 


# Calculate log odds and SE for each group
ORds$logO_BICR<- log(ORds$OR_BICR_exp)
ORds$logO_INV <- log(ORds$OR_INV_exp)

# Compute SE of log odds using 95% CI
ORds$SE_logO_BICR <- (log(ORds$OR_BICR_U_exp) - log(ORds$OR_BICR_L_exp)) / (2 * qnorm(0.975))
ORds$SE_logO_INV <- (log(ORds$OR_INV_U_exp) - log(ORds$OR_INV_L_exp)) / (2 * qnorm(0.975))

# Calculate log OR and its standard error
ORds$logOR <- ORds$logO_INV - ORds$logO_BICR
ORds$SE_logOR <- sqrt(ORds$SE_logO_INV^2 + ORds$SE_logO_BICR^2)

# Perform meta-analysis (random-effects model using REML)
meta<-metagen(TE=logOR,seTE=SE_logOR,data=ORds,sm="OR",common = TRUE,random = FALSE,studlab =PMID)
ss<-summary(meta)
# Convert results to OR scale
pooled_OR <- exp(ss$fixed$TE)
pooled_OR_CI <- exp(c(ss$fixed$lower, ss$fixed$upper))

pooled_OR
pooled_OR_CI

png("forest_plot.png", width = 1000, height = 1200)
forest(meta,backtransf = TRUE)
dev.off()
