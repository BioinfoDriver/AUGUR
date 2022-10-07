
# load data
lihc.risk.score <- readRDS(file='/data/Gao_Cell_2019_risk_score.rds')

# data prepare
cli.sig.char <- lihc.risk.score[, c(1, 4, 3, 10:15, 24, 25, 19:23, 45, 7, 9)]
colnames(cli.sig.char) <- c('patientID', 'age', 'gender', 'cirrhosis', 'tumorNumber', 'tumorSize', 
 'lymphNodeMetastasis', 'tumorThrombus', 'tumourEnapsulation', 'BCLCstage', 'TNMstage', 'totalBilirubin', 
 'ALB', 'ALT', 'GGT', 'AFP', 'risk.categ', 'os_time', 'os')

cli.sig.char$age <- ifelse(cli.sig.char$age >= 60, '≥60', '<60')
cli.sig.char$cirrhosis <- factor(cli.sig.char$cirrhosis, levels=c(0, 1))

cli.sig.char$tumorNumber <- ifelse(cli.sig.char$tumorNumber >= 2, '≥2', '1')
cli.sig.char$tumorNumber <- factor(cli.sig.char$tumorNumber, levels=c('1', '≥2'))

cli.sig.char$tumorSize <- ifelse(cli.sig.char$tumorSize > 5, '>5', '≤5')
cli.sig.char$tumorSize <- factor(cli.sig.char$tumorSize, levels=c('≤5', '>5'))

cli.sig.char$tumourEnapsulation <- factor(cli.sig.char$tumourEnapsulation, levels=c(1, 0))


cli.sig.char$BCLCstage <- ifelse(cli.sig.char$BCLCstage %in% c('C'), 'C', 'A/B')
cli.sig.char$TNMstage <- ifelse(cli.sig.char$TNMstage %in% c('IA', 'IB', 'II'), 'I/II', 'III/IV')
cli.sig.char$totalBilirubin <- ifelse(cli.sig.char$totalBilirubin >= 20.5, 'Yes', 'No')     
cli.sig.char$ALB <- ifelse(cli.sig.char$ALB < 35, 'Yes', 'No')     
cli.sig.char$ALT <- ifelse(cli.sig.char$ALT > 50, 'Yes', 'No')     
cli.sig.char$GGT <- ifelse(cli.sig.char$GGT >= 40, 'Yes', 'No')     
cli.sig.char$AFP <- ifelse(cli.sig.char$AFP > 300, 'Yes', 'No')     

cli.sig.char$risk.categ <- factor(cli.sig.char$risk.categ, levels=c('low risk', 'high risk'))


# 单因素cox分析
UnivariateCox <- function(cli.data, covariates)
{
  library('survival')
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

uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates=colnames(cli.sig.char)[2:17])
# sapply(cli.sig.char[, 2:17], function(x) sum(!is.na(x)))

# 多因素cox分析
cli.sig.char <- cli.sig.char[, c('patientID', 'os', 'os_time', 'risk.categ', 'age', 'gender', 'cirrhosis', 
 'AFP', 'TNMstage', 'BCLCstage', "tumorThrombus", "tumourEnapsulation")]
cli.sig.char <- cli.sig.char[apply(cli.sig.char, 1, function(x) !any(is.na(x))), ]


source('/code/Rscript/Cox.function.R')
uni.mul.cox.res <- Cox.function(time=cli.sig.char$os_time, event=cli.sig.char$os, 
 clinical.data=cli.sig.char, clinical.variate = c(4:9))

