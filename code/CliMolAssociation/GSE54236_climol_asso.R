
# load data
GSE54236.cli.data <- readRDS(file='/data/GSE54236_risk_score.rds')

# HCC patients were separated into more aggressive and less aggressive HCC, according to the tumor doubling time

GSE54236.cli.data$TD <- ifelse(GSE54236.cli.data$doubling_time < 53, 'FastGrowing', 'SlowGrowing')

table(GSE54236.cli.data$TD, GSE54236.cli.data$risk.categ)
fisher.test(table(GSE54236.cli.data$TD, GSE54236.cli.data$risk.categ))$p.value # 0.001146946(max), 0.001146946(aver)
cor.test(GSE54236.cli.data$doubling_time, GSE54236.cli.data$risk.score)$estimate # -0.3559948(max), -0.385938(aver)
cor.test(GSE54236.cli.data$doubling_time, GSE54236.cli.data$risk.score)$p.value # 0.001380078(max), 0.0004835241(aver)


library('ggpubr')

sca.plot <- ggscatter(data=GSE54236.cli.data, x='risk.score', y='doubling_time', color = "black", size = 0.8, cor.coef = TRUE,
 fill = "lightgray", xlab = 'Risk score', ylab = 'HCC doubling time', add = "reg.line")


bar.plot <- ggbarplot(as.data.frame(table(GSE54236.cli.data[, c('TD', 'risk.categ')])), 
 "risk.categ", "Freq", fill = "TD", color="TD", ylab='Number of patients', xlab=FALSE, label=TRUE, 
 lab.col="white", lab.pos="in", lab.size=2, position=position_dodge(0.9),
 title='TD', lab.vjust=2, palette="Paired") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))


ggsave(ggarrange(sca.plot, bar.plot, ncol=4, nrow=4, legend = "none"), 
 file='/result/Section4/GSE54236_max_cli_risk_com.pdf')

