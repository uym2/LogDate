setwd("~/my_gits/LogDate/Data/bio_data/H1N1/")

require(ggplot2)

d = data.frame(method=c("LSD","LF","phyML_B_strict","phyML_B_lnorm","B_strict","B_lnorm","wLogDate"),
               org=c(2008.86,2008.86,2008.19,2007.55,2008.73,2008.50,2008.95),
               med=c(2008.78,2008.82,2008.21,2007.61,2008.74,2008.55,2008.84),
               min=c(2008.13,2008.02,2007.32,2003.54,2008.15,2005.91,2007.81),
               lower=c(2008.34,2008.37,2007.87,2006.77,2008.56,2008.08,2008.18),
               upper=c(2008.95,2008.99,2008.46,2008.16,2008.91,2008.83,2009.05),
               max = c(2009.00,2009.01,2008.63,2008.32,2008.96,2008.93,2009.18))

d$method = factor(d$method,levels=d$method[c(5,6,3,4,1,2,7)])

ggplot(d, aes(x=method,fill=method)) +
  geom_boxplot(
    aes(ymin = min, lower = lower, middle = med, upper = upper, ymax = max),
    stat = "identity" ) + geom_point(aes(y=org)) +
    theme_classic() + theme(legend.position = "bottom",legend.title = element_blank()) +
    theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x=element_blank()) +
    scale_fill_brewer(palette = "Paired")

ggsave("tMRCA.pdf",height = 3,width = 4.5)
ggsave("~/my_gits/LogDate/Figures/Revised_results/Figure_4a_revised.pdf",height = 3,width = 4.5)