# load data
lci.risk.score <- readRDS(file='/data/lci_risk_score.rds')


# data prepare
cli.sig.char <- lci.risk.score[, c("LCS.ID", "Survival.status", "Survival.months", "risk.score", "TNM.staging")]
colnames(cli.sig.char) <- c('patient_id', 'os', 'os_time', 'risk.score', 'pathologic_stage')

cli.sig.char$pathologic_stage[cli.sig.char$pathologic_stage == '.'] <- NA
cli.sig.char$pathologic_stage <- ifelse(cli.sig.char$pathologic_stage %in% c('I', 'II'), 'I/II', 'III/IV')


# nomogram
library(rms)
dd <- datadist(cli.sig.char)
options(datadist='dd')

coxmodel <- cph(Surv(os_time, os==1)~pathologic_stage+risk.score, x=T, y=T, data=cli.sig.char, surv=T)
surv <- Survival(coxmodel)
surv1 <- function(x) surv(1*1*365.25, lp=x)
surv3 <- function(x) surv(1*3*365.25, lp=x)
surv5 <- function(x) surv(1*5*365.25, lp=x)

nom <- nomogram(coxmodel, fun=list(surv1, surv3, surv5), lp=F, 
 funlabel=c('1-Year Survival Prob', '3-Year Survival Prob', '5-Year Survival Prob'), 
 maxscale=100, fun.at=c('1.00', '0.90', '0.70', '0.50', '0.30', '0.10', '0.0'))

setwd('/result/Section6/')
pdf('lci_nom_os.pdf')
 plot(nom)
dev.off()



# Calibration curve 
units(cli.sig.char$os_time) <- 'Day'

oneyear.coxmodel <- cph(Surv(os_time, os==1)~pathologic_stage+risk.score, 
 x=T, y=T, data=cli.sig.char, surv=T, time.inc = 1*365.25) # 5 表示计算1年校准曲线

threeyear.coxmodel <- cph(Surv(os_time, os==1)~pathologic_stage+risk.score, 
 x=T, y=T, data=cli.sig.char, surv=T, time.inc = 3*365.25) # 5 表示计算1年校准曲线

fiveyear.coxmodel <- cph(Surv(os_time, os==1)~pathologic_stage+risk.score, 
 x=T, y=T, data=cli.sig.char, surv=T, time.inc = 5*365.25) # 5 表示计算1年校准曲线


# u 需要与模型中的time.inc一致, m 约等于1/3样本量, 需要不断调试
oneyear.cal <- calibrate(oneyear.coxmodel, cmethod='KM', method='boot', u=1*365.25, m=73, B=1000) 
threeyear.cal <- calibrate(threeyear.coxmodel, cmethod='KM', method='boot', u=3*365.25, m=73, B=1000) 
fiveyear.cal <- calibrate(fiveyear.coxmodel, cmethod='KM', method='boot', u=5*365.25, m=73, B=1000) 



pdf('lci_cali.pdf')
plot(oneyear.cal, lwd=2, lty=1, pch=16, errbar.col='#D9D5D4', xlim=c(0.2, 1), ylim=c(0.2, 1), 
 xlab='Nomogram-predicted Survival', ylab='Proportion of Patients Alive', legend=FALSE,
 col='#D9D5D4', subtitles=F, par.corrected=list(col="#D9D5D4", lty=6, lwd=5, pch=17))

plot(threeyear.cal, lwd=2, lty=1, pch=16, errbar.col='#B0C767', xlim=c(0.2, 1), ylim=c(0.2, 1), 
 xlab='Nomogram-predicted Survival', ylab='Proportion of Patients Alive', legend=FALSE,
 col='#B0C767', subtitles=F, add=TRUE, par.corrected=list(col="#B0C767", lty=1, lwd=2, pch=17))

plot(fiveyear.cal, lwd=2, lty=1, pch=16, errbar.col='#ECD1D2', xlim=c(0.2, 1), ylim=c(0.2, 1), 
 xlab='Nomogram-predicted Survival', ylab='Proportion of Patients Alive', legend=FALSE,
 col='#ECD1D2', subtitles=F, add=TRUE, par.corrected=list(col="#ECD1D2", lty=1, lwd=2, pch=17))

# lines(cal[, c('mean.predicted', 'KM')], type='b', lwd=2, col='red', pch=16)
abline(0, 1, lty=3, lwd=2, col='black')

legend(x = 'bottomright', legend=c('1-year OS', '3-year OS', '5-year OS'), lty=1, lwd=2, 
col = c('#D9D5D4', '#B0C767', '#ECD1D2'), )

dev.off()


