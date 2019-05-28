setwd("/Users/uym2/MyPapers/LogDate/Results")

require(ggplot2)

d = read.table("results_all.txt")
colnames(d) = c("set","rep","method","treeTool","clockModel","benchmark","measure","value")
d$method = factor(d$method,levels=c("BEAST_strict","BEAST_lnorm","lsd","LF","logDate","wlD"))
d$clockModel = factor(d$clockModel,levels=c("lognorm","gamma"))

ggplot(d[d$treeTool != "RAxML" & d$measure=="rmse",],aes(x=clockModel,y=value,group=interaction(method,clockModel),fill=method)) +
  geom_boxplot(outlier.size = 0.5,outlier.alpha = 0.5) + theme_bw() +
  facet_wrap(~benchmark+set,scale="free",nrow = 2) + ylab("rmse") +
  theme(legend.position = "bottom", legend.title = element_blank(),axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")
ggsave("Four_settings.pdf")

ggplot(d[d$treeTool != "RAxML" & d$measure=="rmse" & d$benchmark=="brTime",],aes(x=clockModel,y=value,group=interaction(method,clockModel),fill=method)) +
  #geom_bar(stat="summary",position="dodge",alpha=0.85) + 
  #geom_errorbar(stat="summary",position="dodge",aes(color=method))+ 
  stat_summary(size=0.1,aes(color=method),position=position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~set,scale="free") + ylab("rmse") +
  theme(legend.position = "bottom", legend.title = element_blank(),axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")
ggsave("Four_settings_brTime.pdf")

ggplot(d[d$treeTool != "RAxML" & d$measure=="rmse" & d$benchmark=="nodeAge",],aes(x=clockModel,y=value,group=interaction(method,clockModel),fill=method)) +
  #geom_bar(stat="summary",position="dodge",alpha=0.85) + 
  #geom_errorbar(stat="summary",position="dodge",aes(color=method))+ 
  stat_summary(size=0.1,aes(color=method),position=position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~set,scale="free") + ylab("rmse") +
  theme(legend.position = "bottom", legend.title = element_blank(),axis.title.x = element_blank()) +
  scale_color_brewer(palette = "Dark2")
ggsave("Four_settings_nodeAge.pdf")

ggplot(d[d$treeTool != "RAxML" & d$measure=="rmse",],aes(x=method,y=value,group=method,fill=method)) +
  #geom_bar(stat="summary",alpha=0.85) + 
  #geom_errorbar(stat="summary",width=0.4,aes(color=method)) +
  stat_summary(size=0.1,aes(color=method)) +
  theme_bw() + facet_wrap(~clockModel+benchmark) + ylab("rmse") +
  theme(axis.text.x = element_blank(), legend.position="bottom",axis.title.x=element_blank(),legend.title=element_blank()) + 
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2")
ggsave("Combined.pdf")