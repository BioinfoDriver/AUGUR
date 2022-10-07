
setwd('/data/')
# TCGA
tcga.lihc.risk.score <- readRDS(file='tcga_risk_score.rds') # 2.4959 (1.6953-3.6745) <0.0001
cli.sig.char <- tcga.lihc.risk.score[, c('os', 'os_time', 'risk.categ')]
cli.sig.char$risk.categ <- factor(cli.sig.char$risk.categ, levels=c('low risk', 'high risk'))

# ICGC
icgc.lihc.risk.score <- readRDS(file='icgc_risk_score.rds')# 4.3763 (1.9698-9.7231) 0.000289612650769972
cli.sig.char <- icgc.lihc.risk.score[, c('donor_vital_status', 'donor_survival_time', 'risk.categ')]
colnames(cli.sig.char) <- c('os', 'os_time', 'risk.categ')
cli.sig.char$risk.categ <- factor(cli.sig.char$risk.categ, levels=c('low risk', 'high risk'))

# Gao_Cell_2019
Gao.Cell.2019.risk.score <- readRDS(file='Gao_Cell_2019_risk_score.rds') # 2.4697 (1.408-4.3318) 0.00161353051275268
cli.sig.char <- Gao.Cell.2019.risk.score[, c('Survial  (1, dead; 0, alive)', 'Overall survial (month)', 'risk.categ')]
colnames(cli.sig.char) <- c('os', 'os_time', 'risk.categ')
cli.sig.char$risk.categ <- factor(cli.sig.char$risk.categ, levels=c('low risk', 'high risk'))

# GSE144269
GSE144269.risk.score <- readRDS(file='GSE144269_risk_score.rds') # 2.7426 (1.1269-6.6751) 0.0262045934347918
cli.sig.char <- GSE144269.risk.score[, c('survival.status', 'survival.time', 'risk.categ')]
colnames(cli.sig.char) <- c('os', 'os_time', 'risk.categ')
cli.sig.char$risk.categ <- factor(cli.sig.char$risk.categ, levels=c('low risk', 'high risk'))

# GSE14520
GSE14520.risk.score <- readRDS(file='lci_risk_score.rds') # 2.3305 (1.4912-3.6423) 0.000204106297218611
cli.sig.char <- GSE14520.risk.score[, c('Survival.status', 'Survival.months', 'risk.categ')]
colnames(cli.sig.char) <- c('os', 'os_time', 'risk.categ')
cli.sig.char$risk.categ <- factor(cli.sig.char$risk.categ, levels=c('low risk', 'high risk'))


# GSE1898|4024
lec.risk.score <- readRDS(file='lec_risk_score.rds') # 2.3332 (1.4348-3.794) 0.000636924628483762
cli.sig.char <- lec.risk.score[, c('OS_Status', 'OS_Time', 'risk.categ')]
colnames(cli.sig.char) <- c('os', 'os_time', 'risk.categ')
cli.sig.char$risk.categ <- factor(cli.sig.char$risk.categ, levels=c('low risk', 'high risk'))

# E_TABM_36
E_TABM_36.risk.score <- readRDS(file='E_TABM_36_risk_score.rds') # 3.1788 (1.1999-8.4214) 0.0199896062630961
cli.sig.char <- E_TABM_36.risk.score[, c('os_status', 'os_time', 'risk.categ')]
colnames(cli.sig.char) <- c('os', 'os_time', 'risk.categ')
cli.sig.char$risk.categ <- factor(cli.sig.char$risk.categ, levels=c('low risk', 'high risk'))


   

# 单因素cox分析
UnivariateCox <- function(cli.data, covariates)
{
  #STEP1:构建单因素分析的对象
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(os_time, os)~', x)));
  
  #STEP2:单因素Cox分析
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = cli.data)});

  #STEP3:提取有用信息
  univ_results <- lapply(univ_models, function(x)
                          {                             
                           tmp <-summary(x);
						   
                           #提取p值，保留两位有效数字
                           p.value <- round(tmp$coefficients[ ,5], digits = 4);
                           p.value[which(p.value < 0.0001)] <- "<0.0001";
                           
                           #提取beta值，这里的coefficients为矩阵，但是只有一行，所以可以这样取值
                           #beta <- round(tmp$coefficients[ ,1], digits = 4);
                           
                           #提取风险比
                           HR <- round(tmp$coefficients[ ,2], digits = 4);
                           
                           #提取95%置信区间上下界
                           HR.confint.lower <- round(tmp$conf.int[,"lower .95"], 4);
                           HR.confint.upper <- round(tmp$conf.int[,"upper .95"], 4);    

                           #合并风险比HR和置信区间为一个内容
                           HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")");
                           
                           variate <- rownames(tmp$coefficients);
                           
                           #将所有值合并在一个矩阵中
                           all.data <- as.data.frame(cbind(variate, HR, p.value));
                         }
                        )
  univ_results <- do.call(rbind, univ_results)
  return(univ_results)
}

uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates='risk.categ')



# Function to fit meta-analytic equal-, fixed-, and random-effects models 
# https://stats.stackexchange.com/questions/343316/hazard-ratio-meta-analysis
# https://www.bmj.com/content/343/bmj.d2304
MetaforHrMetaAnalysis <- function(dat, method=c('EE', 'FE', 'REML'), ifplot, fname){
 method <- match.arg(method)
 
 library('metafor')
 library('dplyr')
 
 colnames(dat) <- c('study', 'hr', 'ci.lb', 'ci.ub', 'pval')
 dat <- dat %>% mutate(yi = log(hr), sei = (log(ci.ub)-log(ci.lb))/(2*1.96))
 
 # whether an equal-, a fixed or random-effects model should be fitted.
 res <- rma(yi=yi, sei=sei, data=dat, method=method, slab=study)
 
 meta.sum <- exp(c(res$b[1, 1], res$ci.lb, res$ci.ub))
 meta.sum <- c(meta.sum, res$pval)
 names(meta.sum) <- c('hr', 'ci.lb', 'ci.ub', 'pval')
 
 # plot
 if(ifplot){
  pdf(fname)
   forest(x=res, annotate=TRUE, header=c('Author(s) and Year', 'HR [95% CI]'), refline=1, 
   xlab='Hazard ratio', mlab='Overall', ilab=dat$pval, ilab.xpos=7, ilab.pos=2,
   colout='#3182bd', col='#a50f15', transf=exp)
  dev.off()
 }
 
 return(meta.sum)
}


RmetaHrMetaAnalysis <- function(dat, method=c("fixed", "random"), ifplot, fname){
 method <- match.arg(method)

 library('rmeta')
 library('dplyr')
 
 colnames(dat) <- c('study', 'hr', 'ci.lb', 'ci.ub', 'pval')
 dat <- dat %>% mutate(yi = log(hr), sei = (log(ci.ub)-log(ci.lb))/(2*1.96))
 
 # whether a fixed- or random-effects model should be fitted.
 res <- meta.summaries(d = yi, se = sei, method = method, logscale=TRUE, names=study, data=dat)
 
 meta.sum <- exp(c(res$summary, res$summary - 1.96 * res$se.summary, res$summary + 1.96 * res$se.summary))
 meta.sum <- c(meta.sum, res$test[2])
 names(meta.sum) <- c('hr', 'ci.lb', 'ci.ub', 'pval')
 
 # plot
 if(ifplot){
  pdf(fname)
   metaplot(mn=dat$yi, se=dat$sei, labels=dat$study, xlab='Hazard Ratio', 
   summn = res$summary, sumse = res$se.summary, sumnn= 1/res$se.summary^2, xlim=c(-1, 3), summlabel="Overall",
   zero=0, colors=meta.colors(box="#3182bd",lines="#a50f15", zero="red", summary="black",text="black"), xaxt='n')
   axis(1, at=log(c(0.5,1,2,4,8,16)), labels=c(0.5,1,2,4,8,16))
  dev.off()
 }
 
 return(meta.sum)
}

exam.data <- data.frame(dataset=c('TCGA', 'ICGC', 'Gao et al.', 'GSE144269', 'GSE14520', 'GSE1898|4024', 'E-TABM-36'), 
 HR=c(2.50, 4.38, 2.47, 2.74, 2.33, 2.33, 3.18), 
 lower=c(1.70, 1.97, 1.41, 1.13, 1.49, 1.43, 1.20), upper=c(3.67, 9.72, 4.33, 6.68, 3.64, 3.79, 8.42), 
 pvalue=c('<0.0001', '0.0003', '0.0016', '0.0262', '0.0002', '0.0006', '0.0200'))

RmetaHrMetaAnalysis(exam.data, 'fixed', TRUE, '/result/Section3/metaHRmetafor.pdf')
MetaforHrMetaAnalysis(exam.data, 'FE', TRUE, '/result/Section3/metaHRrmeta.pdf')


