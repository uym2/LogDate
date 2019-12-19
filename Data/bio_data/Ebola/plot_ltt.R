setwd("~/my_gits/LogDate/Data/bio_data/Ebola/")

require("ggplot2")

d = read.table("ltt_data.txt",header = T)

ggplot(d,aes(x=time,y=lineages,color=method)) + geom_line() + theme_classic() +
  scale_color_manual(values=c("#1f78b4","#e31a1c","#fb9a99","#fdbf6f")) +
  theme(legend.position = "None")

ggsave("Ebola_ltt.pdf",height = 3,width = 3)
ggsave("~/my_gits/LogDate/Figures/Revised_results/Figure_4d_revised.pdf",height=3,width=3)