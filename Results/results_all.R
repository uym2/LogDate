setwd("/Users/uym2/MyPapers/LogDate/Results")

require(ggplot2)

d = read.table("results_all.txt")
colnames(d) = c("set","rep","method","treeTool","clockModel","benchmark","measure","value")

ggplot(d[d$treeTool != "RAxML",],aes(x=clockModel,y=value,group=interaction(method,clockModel),fill=method)) +
  geom_boxplot(outlier.size = 0.5,outlier.alpha = 0.5) + theme_bw() +
  facet_wrap(~benchmark+set,scale="free",nrow = 2) + ylab("rmse") +
  theme(legend.position = "bottom", legend.title = element_blank(),axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Accent")

ggsave("Four_settings.pdf")

ggplot(d[d$treeTool != "RAxML",],aes(x=method,y=value,group=method,fill=method)) +
  geom_bar(stat="summary") + geom_errorbar(stat="summary",width=0.4) +
  theme_bw() + facet_wrap(~clockModel+benchmark) + ylab("rmse") +
  theme(axis.text.x = element_blank(), legend.position="bottom",axis.title.x=element_blank(),legend.title=element_blank()) + 
  scale_fill_brewer(palette = "Accent")

ggsave("Combined.pdf")