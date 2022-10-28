library(ggplot2)
library(R.matlab)
library(readr)
library(latex2exp)
library(dplyr)

setwd("~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Privacy_information/Results_epsilon")


k=4
list_epsilon=c(1,10)
alg_string=c("alg4")
# alg_string=c("alg3", "alg4")

df_s1<-c()
for (epsilon in list_epsilon){
  
  for (i in alg_string){
    
    m1 <- read.csv(file = paste("k",as.character(k), i,"_0epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m2 <- read.csv(file = paste("k",as.character(k), i,"_5epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    # m3 <- read.csv(file = paste("k",as.character(k), i,"_3epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    # m4 <- read.csv(file = paste("k",as.character(k), i,"_5epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m5 <- read.csv(file = paste("k",as.character(k), i,"_10epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m6 <- read.csv(file = paste("k",as.character(k), i,"_20epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m7 <- read.csv(file = paste("k",as.character(k), i,"_30epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m8 <- read.csv(file = paste("k",as.character(k), i,"_40epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m9 <- read.csv(file = paste("k",as.character(k), i,"_60epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m10 <- read.csv(file = paste("k",as.character(k), i,"_80epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    
    daf1 <-data.frame(m1)
    daf2 <-data.frame(m2)
    # daf3 <-data.frame(m3)
    # daf4 <-data.frame(m4)
    daf5 <-data.frame(m5)
    daf6 <-data.frame(m6)
    daf7 <-data.frame(m7)
    daf8 <-data.frame(m8)
    daf9 <-data.frame(m9)
    daf10 <-data.frame(m10)
    
    
    
    
    colnames(daf1) <- c('t')
    colnames(daf2) <- c('t')
    # colnames(daf3) <- c('t')
    # colnames(daf4) <- c('t')
    colnames(daf5) <- c('t')
    colnames(daf6) <- c('t')
    colnames(daf7) <- c('t')
    colnames(daf8) <- c('t')
    colnames(daf9) <- c('t')
    colnames(daf10) <- c('t')
    
    meants1= mean(daf1$t)
    meants2=mean(daf2$t)
    # meants3=mean(daf3$t)
    # meants4=mean(daf4$t)
    meants5=mean(daf5$t)
    meants6= mean(daf6$t)
    meants7=mean(daf7$t)
    meants8=mean(daf8$t)
    meants9=mean(daf9$t)
    meants10=mean(daf10$t)
    
    sdts1=sd(daf1$t)
    sdts2=sd(daf2$t)
    # sdts3=sd(daf3$t)
    # sdts4=sd(daf4$t)
    sdts5=sd(daf5$t)
    sdts6=sd(daf6$t)
    sdts7=sd(daf7$t)
    sdts8=sd(daf8$t)
    sdts9=sd(daf9$t)
    sdts10=sd(daf10$t)
    
    
    
    Time <- c(meants1,meants2, meants5,meants6,meants7,meants8,meants9, meants10)
    std_dev <- c(sdts1,sdts2,sdts5,sdts6,sdts7,sdts8,sdts9,sdts10)
    algo= rep(i,8)
    eps=rep(epsilon,8)
    L<-c(0,5,10,20,30,40,60,80)
    
    
    df4 <- data.frame(Time, std_dev, algo,L,eps)
    df4$algo<-as.factor(df4$algo)
    df4$eps<-as.factor(df4$eps)
    df_s1<-rbind(df_s1,df4)
  }
  
}

list_epsilon=c(0.1,1)
alg_string=c("alg3")


df_s2<-c()
for (epsilon in list_epsilon){
  
  for (i in alg_string){
    
    m1 <- read.csv(file = paste("k",as.character(k), i,"_0epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m2 <- read.csv(file = paste("k",as.character(k), i,"_5epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    # m3 <- read.csv(file = paste("k",as.character(k), i,"_3epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    # m4 <- read.csv(file = paste("k",as.character(k), i,"_5epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m5 <- read.csv(file = paste("k",as.character(k), i,"_10epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m6 <- read.csv(file = paste("k",as.character(k), i,"_20epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m7 <- read.csv(file = paste("k",as.character(k), i,"_30epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m8 <- read.csv(file = paste("k",as.character(k), i,"_40epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m9 <- read.csv(file = paste("k",as.character(k), i,"_60epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m10 <- read.csv(file = paste("k",as.character(k), i,"_80epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    
    daf1 <-data.frame(m1)
    daf2 <-data.frame(m2)
    # daf3 <-data.frame(m3)
    # daf4 <-data.frame(m4)
    daf5 <-data.frame(m5)
    daf6 <-data.frame(m6)
    daf7 <-data.frame(m7)
    daf8 <-data.frame(m8)
    daf9 <-data.frame(m9)
    daf10 <-data.frame(m10)
    
    
    
    
    colnames(daf1) <- c('t')
    colnames(daf2) <- c('t')
    # colnames(daf3) <- c('t')
    # colnames(daf4) <- c('t')
    colnames(daf5) <- c('t')
    colnames(daf6) <- c('t')
    colnames(daf7) <- c('t')
    colnames(daf8) <- c('t')
    colnames(daf9) <- c('t')
    colnames(daf10) <- c('t')
    
    meants1= mean(daf1$t)
    meants2=mean(daf2$t)
    # meants3=mean(daf3$t)
    # meants4=mean(daf4$t)
    meants5=mean(daf5$t)
    meants6= mean(daf6$t)
    meants7=mean(daf7$t)
    meants8=mean(daf8$t)
    meants9=mean(daf9$t)
    meants10=mean(daf10$t)
    
    sdts1=sd(daf1$t)
    sdts2=sd(daf2$t)
    # sdts3=sd(daf3$t)
    # sdts4=sd(daf4$t)
    sdts5=sd(daf5$t)
    sdts6=sd(daf6$t)
    sdts7=sd(daf7$t)
    sdts8=sd(daf8$t)
    sdts9=sd(daf9$t)
    sdts10=sd(daf10$t)
    
    
    
    Time <- c(meants1,meants2, meants5,meants6,meants7,meants8,meants9, meants10)
    std_dev <- c(sdts1,sdts2,sdts5,sdts6,sdts7,sdts8,sdts9,sdts10)
    algo= rep(i,8)
    eps=rep(epsilon,8)
    L<-c(0,5,10,20,30,40,60,80)
    
    
    df4 <- data.frame(Time, std_dev, algo,L,eps)
    df4$algo<-as.factor(df4$algo)
    df4$eps<-as.factor(df4$eps)
    df_s1<-rbind(df_s1,df4)
  }
  
}


######################################################


df<-rbind(df_s1, df_s2)

colnames(df) <- c('Time', 'std_dev', 'Algorithm', 'L', 'Epsilon')


# df<- df %>% 
#   mutate(Algorithm = ifelse(as.character(Algorithm) == "alg0", "Greedy Alg.", as.character(Algorithm)))

df<-df %>% 
  mutate(Algorithm = ifelse(as.character(Algorithm) == "alg3", "Exponential Mechanism", as.character(Algorithm)))

df<-df %>% 
  mutate(Algorithm = ifelse(as.character(Algorithm) == "alg4", "Randomized Response", as.character(Algorithm)))



g<- ggplot(df, aes(x=L, y=Time, colour=Algorithm, fill=Algorithm)) + 
  geom_errorbar(aes(ymin=Time-(1.96/(sqrt(800)))*std_dev, ymax=Time+(1.96/(sqrt(800)))*std_dev), width=.1) +
  geom_line(aes(group=interaction(Epsilon, Algorithm))) +
  geom_point(aes(shape=(Epsilon), size=0.5), size=2)+
  geom_hline(aes(yintercept=422, linetype="Greedy Algorithm"),color ="#C0392B")+
  theme_light()+
  guides(shape = guide_legend(title = expression(paste(epsilon))))+
  scale_color_manual(values=c("#F59E07","#1F618D", "#1F618D"))+
  labs(x = "Number of influece samples", y="Expected size of the spread")+
  scale_linetype_manual(name = "", values = c(1), 
                        guide = guide_legend(override.aes = list(color = c("#C0392B"))))+
  theme()

g


saving_name="All_comparison.png"

ggsave(saving_name)

###================= code to plot With vs Without post processing



k=4
list_epsilon=c(1,10)
alg_string=c("alg4", "alg6")
# alg_string=c("alg3", "alg4")

df<-c()
for (epsilon in list_epsilon){
  
  for (i in alg_string){
    
    m1 <- read.csv(file = paste("k",as.character(k), i,"_0epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m2 <- read.csv(file = paste("k",as.character(k), i,"_5epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    # m3 <- read.csv(file = paste("k",as.character(k), i,"_3epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    # m4 <- read.csv(file = paste("k",as.character(k), i,"_5epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m5 <- read.csv(file = paste("k",as.character(k), i,"_10epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m6 <- read.csv(file = paste("k",as.character(k), i,"_20epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m7 <- read.csv(file = paste("k",as.character(k), i,"_30epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m8 <- read.csv(file = paste("k",as.character(k), i,"_40epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m9 <- read.csv(file = paste("k",as.character(k), i,"_60epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    m10 <- read.csv(file = paste("k",as.character(k), i,"_80epsilon_",as.character(epsilon),".csv",sep=""), header=F)
    
    daf1 <-data.frame(m1)
    daf2 <-data.frame(m2)
    # daf3 <-data.frame(m3)
    # daf4 <-data.frame(m4)
    daf5 <-data.frame(m5)
    daf6 <-data.frame(m6)
    daf7 <-data.frame(m7)
    daf8 <-data.frame(m8)
    daf9 <-data.frame(m9)
    daf10 <-data.frame(m10)
    
    
    
    
    colnames(daf1) <- c('t')
    colnames(daf2) <- c('t')
    # colnames(daf3) <- c('t')
    # colnames(daf4) <- c('t')
    colnames(daf5) <- c('t')
    colnames(daf6) <- c('t')
    colnames(daf7) <- c('t')
    colnames(daf8) <- c('t')
    colnames(daf9) <- c('t')
    colnames(daf10) <- c('t')
    
    meants1= mean(daf1$t)
    meants2=mean(daf2$t)
    # meants3=mean(daf3$t)
    # meants4=mean(daf4$t)
    meants5=mean(daf5$t)
    meants6= mean(daf6$t)
    meants7=mean(daf7$t)
    meants8=mean(daf8$t)
    meants9=mean(daf9$t)
    meants10=mean(daf10$t)
    
    sdts1=sd(daf1$t)
    sdts2=sd(daf2$t)
    # sdts3=sd(daf3$t)
    # sdts4=sd(daf4$t)
    sdts5=sd(daf5$t)
    sdts6=sd(daf6$t)
    sdts7=sd(daf7$t)
    sdts8=sd(daf8$t)
    sdts9=sd(daf9$t)
    sdts10=sd(daf10$t)
    
    
    
    Time <- c(meants1,meants2, meants5,meants6,meants7,meants8,meants9, meants10)
    std_dev <- c(sdts1,sdts2,sdts5,sdts6,sdts7,sdts8,sdts9,sdts10)
    algo= rep(i,8)
    eps=rep(epsilon,8)
    L<-c(0,5,10,20,30,40,60,80)
    
    
    df4 <- data.frame(Time, std_dev, algo,L,eps)
    df4$algo<-as.factor(df4$algo)
    df4$eps<-as.factor(df4$eps)
    df<-rbind(df,df4)
  }
  
}



######################################################




colnames(df) <- c('Time', 'std_dev', 'Algorithm', 'L', 'Epsilon')


# df<- df %>% 
#   mutate(Algorithm = ifelse(as.character(Algorithm) == "alg0", "Greedy Alg.", as.character(Algorithm)))

df<-df %>% 
  mutate(Algorithm = ifelse(as.character(Algorithm) == "alg6", "RR without post proc.", as.character(Algorithm)))

df<-df %>% 
  mutate(Algorithm = ifelse(as.character(Algorithm) == "alg4", "RR with post proc.", as.character(Algorithm)))



g<- ggplot(df, aes(x=L, y=Time, colour=Algorithm, fill=Algorithm)) + 
  geom_errorbar(aes(ymin=Time-(1.96/(sqrt(800)))*std_dev, ymax=Time+(1.96/(sqrt(800)))*std_dev), width=.1) +
  geom_line(aes(group=interaction(Epsilon, Algorithm))) +
  geom_point(aes(shape=(Epsilon), size=0.5), size=2)+
  geom_hline(aes(yintercept=422, linetype="Greedy Algorithm"),color ="#C0392B")+
  theme_light()+
  guides(shape = guide_legend(title = expression(paste(epsilon))))+
  scale_color_manual(values=c("#F59E07","#1F618D", "#1F618D"))+
  labs(x = "Number of influece samples", y="Expected size of the spread")+
  scale_linetype_manual(name = "", values = c(1), 
                        guide = guide_legend(override.aes = list(color = c("#C0392B"))))+
  theme()

g


saving_name="All_comparison.png"

ggsave(saving_name)



