
# load data
lci.risk.score <- readRDS(file='/data/lci_risk_score.rds')

lci.risk.score <- lci.risk.score[, c(3, 5, 	8:18, 24, 19, 20)]
colnames(lci.risk.score) <- c('Affy_GSM',  'Metastasis_risk', 'Gender', 'Age', 'HBV_status', 'ALT', 'Tumor_size', 
 'Multinodular', 'Cirrhosis', 'TNM_staging', 'BCLC_staging', 'CLIP_staging', 'AFP', 'risk.categ', 'os', 'os_time')

# data prepare
lci.risk.score$Metastasis_risk <- factor(lci.risk.score$Metastasis_risk, levels=c('low', 'high'))
lci.risk.score$Age <- ifelse(lci.risk.score$Age < 60, '<60', '≥60')
lci.risk.score$ALT <- factor(lci.risk.score$ALT, levels=c('low', 'high'))

lci.risk.score$Tumor_size[lci.risk.score$Tumor_size == '.'] <- NA
lci.risk.score$Tumor_size <- factor(lci.risk.score$Tumor_size, levels=c('small', 'large'))

lci.risk.score$TNM_staging[lci.risk.score$TNM_staging == '.'] <- NA
lci.risk.score$TNM_staging <- ifelse(lci.risk.score$TNM_staging %in% c('I', 'II'), 'I/II', 'III/IV')
lci.risk.score$BCLC_staging[lci.risk.score$BCLC_staging == '.'] <- NA
lci.risk.score$BCLC_staging <- ifelse(lci.risk.score$BCLC_staging %in% c('0', 'A'), '0/A', 'B/C')
lci.risk.score$CLIP_staging[lci.risk.score$CLIP_staging == '.'] <- NA
lci.risk.score$CLIP_staging <- ifelse(lci.risk.score$CLIP_staging %in% c('0', '1'), '0/1', '≥2')
lci.risk.score$AFP[lci.risk.score$AFP == '.'] <- NA

lci.risk.score$AFP <- factor(lci.risk.score$AFP, levels=c('low', 'high'))
lci.risk.score$risk.categ <- factor(lci.risk.score$risk.categ, levels=c('low risk', 'high risk'))


cli.sig.char <- lci.risk.score


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

uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates=colnames(cli.sig.char)[2:14])
# sapply(cli.sig.char[, 2:14], function(x) sum(!is.na(x)))

# 多因素cox分析
cli.sig.char <- cli.sig.char[, c('Affy_GSM', 'os', 'os_time', 'risk.categ', 'Age', 'Gender', 
 'Cirrhosis', 'AFP', "TNM_staging", "BCLC_staging", "CLIP_staging")]
cli.sig.char <- cli.sig.char[apply(cli.sig.char, 1, function(x) !any(is.na(x))), ]


source('/code/Rscript/Cox.function.R')
uni.mul.cox.res <- Cox.function(time=cli.sig.char$os_time, event=cli.sig.char$os, 
 clinical.data=cli.sig.char, clinical.variate = c(4:9))
