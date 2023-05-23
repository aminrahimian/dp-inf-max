# Plot figure 2 

library(ggplot2)
library(R.matlab)
library(readr)
library(latex2exp)
library(dplyr)
library(tidyr)


df1<-read.csv('exp_mech.csv')

df2<-read.csv('randomized_version.csv')

df3<- read.csv('randomized_version_wpp.csv')


df1 <-df1 %>%mutate(algorithm= rep('exponential\nmechanism',dim(.)[1])) %>%  
  mutate_at(c('m', 'k','epsilon','algorithm'), as.factor)

df2<- df2 %>%mutate(algorithm= rep('randomize\nresponse',dim(.)[1])) %>%  
  mutate_at(c('m', 'k','epsilon','algorithm'), as.factor)

df3<- df3 %>%mutate(algorithm= rep('rr w/o\npost. proces.',dim(.)[1])) %>%  
  mutate_at(c('m', 'k','epsilon','algorithm'), as.factor) 


df<-rbind(df1,df2,df2,df3)

df <- df %>% mutate(comp=rep(c('c1','c2'),each=dim(df2)[1]*2))
df$comp<- as.factor(df$comp)


df<- df %>% mutate(s=paste('k = ',k, sep = ""))

df$s<- factor(df$s, levels = c('k = 4', 'k = 8', 'k = 12'), ordered = TRUE) 


str(df)

g<- ggplot(df, aes(x=m, y=ixs_mu, colour=algorithm, fill=algorithm)) + 
  geom_errorbar(aes(ymin=ixs_mu-(1.96/(sqrt(ixs_n)))*ixs_sd, ymax=ixs_mu+(1.96/(sqrt(ixs_n)))*ixs_sd), width=.05) +
  geom_line(aes(group=interaction(epsilon, algorithm))) +
  geom_point(aes(shape=(epsilon), size=0.5), size=2)+
  theme_light() +
  theme(legend.key.height=unit(0.8, "cm"), legend.key.width =unit(0.1, "cm"),
        legend.text=element_text(size=10))+
  scale_color_manual(values=c("#F59E07","#1F618D", "#C0392B"))+
  labs(x = "number of influece samples", y="expected size of the spread")+
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.x = element_text(colour = "black",size = 10),
  )+
  facet_grid(comp~s)+
  guides(shape = guide_legend(title = expression(paste(epsilon))))

g
