
HmiscCIndex <- function(clinical.data, time, event, models, diff = TRUE){
  
  options(stringsAsFactors = FALSE)
  
  # load survival, Hmisc package
  suppressPackageStartupMessages(require(survival))
  suppressPackageStartupMessages(require(Hmisc))
  
  ##-------------------------------
  # CoxphObject：create a coxph object
  CoxphObject <- function(data, time, event, variables){
    formula <- as.formula(paste("Surv(time, event)~", paste(variables, collapse=" + ")))
    coxph.object <-coxph(formula, data = data)
    
    return(coxph.object)
  }
  
  ##-------------------------------
  # 输出的Cindex和置信区间
  c.ci <- function(rcorrobj){
    CIndex <- sprintf('%0.2f', rcorrobj['C Index'])
    se     <- rcorrobj['S.D.']/2
    Lower <- sprintf('%0.2f', rcorrobj['C Index'] - 1.96*se)
    Upper <- sprintf('%0.2f', rcorrobj['C Index'] + 1.96*se)

    result <- c(CIndex, Lower, Upper)
    names(result) <- c("C-Index", "Lower", "Upper")
    
    return(result)
  }
  
  # 计算每个模型的C-index值和置信区间
  coxph.list <- lapply(models, function(x){CoxphObject(data = clinical.data, time = time, event = event, variables = x)})
  pred.models.coxph <- lapply(coxph.list, function(x){predict(x, type = "risk")}) 
  models.result <- lapply(pred.models.coxph, function(x){rcorr.cens(-x, Surv(time = time, event = event))})
  
  models.filter.result <- lapply(models.result, function(x){c.ci(rcorrobj = x)})
  
  # 是否进行C-index1的比较
  if (diff == FALSE) {
    result <- do.call(rbind, models.filter.result)
    conf.interval <- paste(result[, 'Lower'], result[, 'Upper'], sep = "-")
    result <- data.frame(result[, 'C-Index'], conf.interval)
    colnames(result) <- c('C-index', '95%CI')
    
    return(result)
  } else {
    # 计算比较的p值，都是和第一个模型比较
    compare.cindex <- lapply(pred.models.coxph[-1], function(x){rcorrp.cens(pred.models.coxph[[1]], x, Surv(time = time, event = event))})
    p.value <- c("-", unlist(lapply(compare.cindex, function(x)(round(2*(1 - pnorm(abs(x[['Dxy']] / x[['S.D.']]))), digits=4))))) 

	
    result <- do.call(rbind, models.filter.result)
    conf.interval <- paste(result[, 'Lower'], result[, 'Upper'], sep = "-")
    result <- data.frame(result[, 'C-Index'], conf.interval, p.value)
    colnames(result) <- c('C-index', '95%CI', 'p-value')
    
    return(result)
  }
}  


#  calculate the C-index for a cox model
load(file='/data/four_dataset_integ_model_socre.RData')
risk.sig <- list(c('risk.score', 'pathologic_stage_category'), 'risk.score', 'pathologic_stage_category')
HmiscCIndex(tcga.model.score, tcga.model.score$os_time, tcga.model.score$os, risk.sig, diff = TRUE)

risk.sig <- list(c('risk.score', 'pathologic_stage'), 'risk.score', 'pathologic_stage')
HmiscCIndex(icgc.model.score, icgc.model.score$os_time, icgc.model.score$os, risk.sig, diff = TRUE)

HmiscCIndex(gao.etal.model.score, gao.etal.model.score$os_time, gao.etal.model.score$os, risk.sig, diff = TRUE)

HmiscCIndex(lci.model.score, lci.model.score$os_time, lci.model.score$os, risk.sig, diff = TRUE)
