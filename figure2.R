# Plot expected spread size for dp algorithms for a given size of seed set k
#setwd("~/PycharmProjects/dp-inf-max")

library(ggplot2)
library(R.matlab)
library(readr)
library(latex2exp)
library(dplyr)
library(tidyr)


df1<-read.csv('./data_msm_network/exp_mech.csv')
df2<-read.csv('./data_msm_network/randomized_version.csv')
df3 <- read.csv('./data_msm_network/greedy_alg.csv')
greedy_vals=read.csv('./data_msm_network/greedy_alg_reference.csv')


df1 <-df1 %>%mutate(algorithm= rep('exponential\nmechanism',dim(.)[1])) %>%  
  mutate_at(c('m', 'k','epsilon','algorithm'), as.factor) 


df2<- df2 %>%mutate(algorithm= rep('randomized\nresponse',dim(.)[1])) %>%  
  mutate_at(c('m', 'k','epsilon','algorithm'), as.factor) 
 


df3<- df3 %>%mutate(algorithm= rep('non-private\nalgorithm',dim(.)[1])) %>%  
  mutate_at(c('m', 'k','epsilon','algorithm'), as.factor) %>% 
  mutate(epsilon=NA) 


df<-rbind(df1,df2,df3)

df<- df %>% mutate(s=paste('k = ',k, sep = ""))
df$s<- factor(df$s, levels = c('k = 10', 'k = 20'), ordered = TRUE) 
hline_dat = data.frame(s=c('k = 10', 'k = 20'),
                       threshold=c(greedy_vals[1,2],greedy_vals[2,2]))
df<-df %>% mutate(m=as.numeric(as.character(df$m)))


g<- ggplot(df, aes(x=m, y=ixs_mu, colour=algorithm)) + 
  geom_errorbar(aes(ymin=ixs_mu-(1.96/(sqrt(ixs_n)))*ixs_sd, ymax=ixs_mu+(1.96/(sqrt(ixs_n)))*ixs_sd), width=.05) +
  geom_point(aes(shape=(epsilon), size=0.5), size=2,na.rm =T)+
  geom_line(aes(group=interaction(epsilon, algorithm))) +
  geom_hline(data=hline_dat, aes(yintercept=threshold,linetype = "greedy baseline"), 
             colour="#A21212")+
  scale_x_continuous(limits = c(0, 3000)) +
  theme_light()+
  scale_color_manual(values=c("#F59E07","#1F618D", "#C0392B"))+
  labs(x = "number of influece samples", y="expected spread  size")+
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.x = element_text(colour = "black",size = 10),
  )+
  facet_grid(.~s)+
  guides(shape = guide_legend(title = expression(paste(epsilon))))+
  scale_shape(na.translate = FALSE)+
  scale_linetype_manual(name = "", values = c(2), 
                        guide = guide_legend(override.aes = list(color = c("#A21212"))))+
  theme(legend.key.height=unit(0.8, "cm"), legend.key.width =unit(0.1, "cm"),
        legend.text=element_text(size=9))
g



saving_name="HIV-network-k-1020-exp-rr.pdf"


ggsave(saving_name, plot = g, width = 6.5, height = 4, units = "in", dpi = 300)


