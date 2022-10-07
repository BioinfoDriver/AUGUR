
# load data
hcc.cli.data <- readRDS(file='/data/E_TABM_36_risk_score.rds')


hcc.cli.data$Age <- ifelse(hcc.cli.data$Age < 60, '<60', '≥60')
hcc.cli.data$Age <- factor(hcc.cli.data$Age, levels=c('<60', '≥60'))
hcc.cli.data$HBV_status <- ifelse(hcc.cli.data$HBV_status == 'HBV_titer_negative', 'No', 'Yes')
hcc.cli.data$risk.categ <- factor(hcc.cli.data$risk.categ, levels=c('low risk', 'high risk'))
hcc.cli.data$DiseaseStage[hcc.cli.data$DiseaseStage == '  '] <- NA
colnames(hcc.cli.data)[11:12] <- c('os', 'os_time')


cli.sig.char <- hcc.cli.data


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

uni.cox.res <- UnivariateCox(cli.data = cli.sig.char, covariates=colnames(cli.sig.char)[c(3:5, 8, 16)])
# sapply(cli.sig.char[, c(3:5, 8, 16)], function(x) sum(!is.na(x)))

# 多因素cox分析
cli.sig.char <- cli.sig.char[, c('PatientID', 'os', 'os_time', 'risk.categ', 'Age', 'Sex')] # , 'DiseaseStage', 'HBV_status'
cli.sig.char <- cli.sig.char[apply(cli.sig.char, 1, function(x) !any(is.na(x))), ]



source('/code/Rscript/Cox.function.R')
uni.mul.cox.res <- Cox.function(time=cli.sig.char$os_time, event=cli.sig.char$os, 
 clinical.data=cli.sig.char, clinical.variate = c(4:6))

