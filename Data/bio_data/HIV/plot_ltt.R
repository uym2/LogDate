setwd("~/my_gits/LogDate/Data/bio_data/HIV/")

require("ggplot2")

d = read.table("ltt_data.txt",header = T)

ggplot(d,aes(x=time,y=lineages,color=method)) + geom_line() + theme_classic() +
  scale_color_manual(values=c("#e31a1c","#fb9a99","#fdbf6f")) +
  theme(legend.position = "None")

ggsave("HIV_ltt.pdf",height = 3,width = 3)