setwd("~/my_gits/LogDate/Data/bio_data/H1N1/")

require("ggplot2")

d = read.table("ltt_data.txt",header = T)

d$method = factor(d$method,levels=levels(d$method)[c(2,1,6,5,4,3,7)])

ggplot(d,aes(x=time,y=lineages,color=method)) + geom_line() + theme_classic() +
  theme(legend.position = "None") + 
  scale_color_brewer(palette = "Paired")

ggsave("H1N1_ltt.pdf",height = 3,width = 3)
ggsave("~/my_gits/LogDate/Figures/Revised_results/Figure_4b_revised.pdf",height=3,width=3)