library(GEOquery)
workdir = "C:/Users/thinkpad/Desktop/design/data/expresstiondata/GSE52562"
setwd(workdir)

gse52562 <- getGEO(filename='GSE52562_series_matrix(1).txt',GSEMatrix =TRUE,getGPL=F)
gpl10558<-getGEO('GPL10558',destdir = ".")
exprs(gse52562)
colnames<-names(pData(phenoData(gse52562)))

library(genefilter)
#func<-ttest(match52562$Symbol,p=0.1)
#gse52562_filter<-genefilter(match52562,filterfun(func)) 
gse52562.filter <- nsFilter(gse52562, require.entrez=F, remove.dupEntrez=F)
gse52562.filter$filter.log

a<-pData(phenoData(gse52562))
b<-a[which(a$`treatment:ch1`=='pre-pidilizumab'),]
filter1<-rownames(b)
pheno52562<-b[,c("pfs.days:ch1",'pfs.status.censorship:ch1',"age:ch1")]
hpd<-pheno52562[c('GSM1269873 ','GSM1269893'),]
non_hpd<-pheno52562[-which(rownames(pheno52562) %in% c('GSM1269873','GSM1269893')),]

exprs52562_f<-exprs(gse52562.filter$eset)[,filter1]
exprs52562_f<-as.data.frame(exprs52562_f)
genenames<-Table(gpl10558)[,c("ID","Symbol")]

exprs52562_f$ID <- rownames(exprs52562_f)
match52562_f<-merge(x=exprs52562_f,y=genenames,by='ID',all.x=T)
match52562_f<-match52562_f[,-1]

failed<-match52562_f[match52562_f$Symbol=="",]
match52562_f<-match52562_f[match52562_f$Symbol!='',]

match52562_f_sum <- aggregate(x = match52562_f,by = list(match52562_f$Symbol), FUN = max)#保留最大值
match52562_f_sum <-match52562_f_sum[,-1]

#COXPH
library(survival)
library(survminer)
test<-cbind(t(match52562_f_sum[,-19]),pheno52562)
colnames(test)[1:17371]<-match52562_f_sum$Symbol
colnames(test)[17372:17374]<-c('time','status','age')
colnames(test)[1:8]<-c('MAR1','MAR2','MAR3','MAR5','MAR6','MAR7','MAR8','MAR9')#处理一下被命名错的基因
test[,'time']<-as.numeric(test[,'time'])
test[,'status']<-as.numeric(test[,'status'])
test.cox <- coxph(Surv(time,status) ~ BANP, data = test)
#同时处理多个基因
covariates <- as.character(colnames(test)[1:17371]) 
covariates_1 <- gsub('-','',covariates)
test_1<-test
colnames(test_1)<-gsub('-','',colnames(test))
colnames(test_1)<-gsub(' ','',colnames(test_1))
univ_formulas <- sapply(covariates_1,function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = test_1)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         z<-x$coefficients[4]
                         res<-c(beta, HR, wald.test, p.value,z)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value","z")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
as.data.frame(res[,which(res$p.value<0.05)])
#not necessary
res_1[,2]<-as.numeric(res[,3])
res_1[,3]<-as.numeric(res[,4])
res_1[,4]<-as.numeric(res[,5])
res_1<-res_1[,-5]
colnames(res_1)<-c('beta','wald.test','p.value','z')
u<-dflme1[,2:60]<-lapply(dflme1[,2:60],as.numeric)
