
#  calculate the C-index for a cox model
#'@description 计算生存模型的C-index或者两个模型比较的p值
#'@param clinical.data: 临床变量信息，至少包含model中的变量信息
#'@param time: 数值型向量，患者对应的生存时间
#'@param event: 患者对应的生存状态，通常0=alive, 1=dead
#'@param models：列表，每一个元素是一个字符型向量，包含一个model所有变量对应的列名,至少包含两个元素
#'@param diff: 逻辑值，diff = FALSE不计算两个模型比较的p值，只返回对应的C-index值；diff = TRUE返回模型对应的C-index值及模型比较的p值
#'@return 返回一个data.frame，包含C-index值，置信区间，比较的p值（当diff = T时）
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


SurvcompCIndex <- function(clinical.data, time, event, models, diff = TRUE){
  
  options(stringsAsFactors = FALSE)
  
  # load survival, survcomp package
  suppressPackageStartupMessages(require(survival))
  suppressPackageStartupMessages(require(survcomp))
  
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
    CIndex <- sprintf('%0.2f', rcorrobj$c.index)
    Lower <- sprintf('%0.2f', rcorrobj$lower)
    Upper <- sprintf('%0.2f', rcorrobj$upper)
    result <- c(CIndex, Lower, Upper)
    names(result) <- c("CIndex", "Lower", "Upper")
    
    return(result)
  }
  
  # 计算每个模型的C-index值和置信区间
  coxph.list <- lapply(models, function(x){CoxphObject(data = clinical.data, time = time, event = event, variables = x)})
  pred.models.coxph <- lapply(coxph.list, function(x){predict(x, type = "risk")}) 
  models.result <- lapply(pred.models.coxph, function(x){concordance.index(x, surv.time= time, surv.event = event, method = "noether")})
  
  models.filter.result <- lapply(models.result, function(x){c.ci(rcorrobj = x)})
  
  # 是否进行C-index1的比较
  if (diff == FALSE) {
    result <- as.data.frame(do.call(rbind, models.filter.result))
    conf.interval <- paste0(result$CIndex, '(', result$Lower, '-', result$Upper, ')')
    result <- data.frame(result$CIndex, conf.interval)
    colnames(result) <- c('CIndex', '95%CI')
    
    return(result)
  } else {
    # 计算比较的p值，都是和第一个模型比较
    compare.cindex.p.value <- sapply(models.result[-1], function(x){cindex.comp(models.result[[1]], x)$p.value})
    
	p.value <- c("-", round(compare.cindex.p.value, digits=4)) 

	
    result <- as.data.frame(do.call(rbind, models.filter.result))
    conf.interval <- paste0(result$CIndex, '(', result$Lower, '-', result$Upper, ')')
    result <- data.frame(result$CIndex, conf.interval, p.value)
    colnames(result) <- c('CIndex', '95%CI', 'p-value')
    
    return(result)
  }
}  



clonal.risk.score <- readRDS(file='/data/tcga_risk_score.rds')
signature.risk.score <- readRDS(file='/data/curated_sig_tcga_risk_score.rds')
colnames(signature.risk.score) <- gsub('\\:', '_', colnames(signature.risk.score))

rownames(signature.risk.score) <- substr(rownames(signature.risk.score), 1, 12)
clinical.data <- cbind(clonal.risk.score, signature.risk.score[clonal.risk.score$patient_id, ])

# saveRDS(clinical.data, file='/data/tcga_sig_risk_score.rds')

risk.sig <- c('risk.score', 'PMID_32198063', 'PMID_30885723', 'PMID_35123387', 'PMID_31335995', 'PMID_31695766',
 'PMID_33033585', 'PMID_34676211', 'PMID_34277622', 'PMID_33777771', 'PMID_34975331', 'PMID_25666192','PMID_22105560')

hmisc.cindex <- HmiscCIndex(clinical.data, clinical.data$os_time, clinical.data$os, risk.sig, diff = TRUE)

survcomp.cindex <- SurvcompCIndex(clinical.data, clinical.data$os_time, clinical.data$os, risk.sig, diff = TRUE)


