library(ggplot2)
library(R.matlab)
library(readr)
library(latex2exp)

setwd("~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Privacy_information/Results")


k=1
epsilon=0.1
alg_string=c("alg3", "alg4")


df<-c()

for (i in alg_string){

  m2_eps01 <- read.csv(file = paste("k",as.character(k), i,"_200epsilon_",as.character(epsilon),".csv",sep=""), header=F)
  m5_eps01 <- read.csv(file = paste("k",as.character(k), i,"_500epsilon_",as.character(epsilon),".csv",sep=""), header=F)
  m10_eps01 <- read.csv(file = paste("k",as.character(k), i,"_1000epsilon_",as.character(epsilon),".csv",sep=""), header=F)
  m15_eps01 <- read.csv(file = paste("k",as.character(k), i,"_1500epsilon_",as.character(epsilon),".csv",sep=""), header=F)
  m20_eps01 <- read.csv(file = paste("k",as.character(k), i,"_2000epsilon_",as.character(epsilon),".csv",sep=""), header=F)
  
  
  daf1 <-data.frame(m2_eps01)
  daf2 <-data.frame(m5_eps01)
  daf3 <-data.frame(m10_eps01)
  daf4 <-data.frame(m15_eps01)
  daf5 <-data.frame(m20_eps01)
  
  
  colnames(daf1) <- c('t')
  colnames(daf2) <- c('t')
  colnames(daf3) <- c('t')
  colnames(daf4) <- c('t')
  colnames(daf5) <- c('t')
  
  meants1= mean(daf1$t)
  meants2=mean(daf2$t)
  meants3=mean(daf3$t)
  meants4=mean(daf4$t)
  meants5=mean(daf5$t)
  
  sdts1=sd(daf1$t)
  sdts2=sd(daf2$t)
  sdts3=sd(daf3$t)
  sdts4=sd(daf4$t)
  sdts5=sd(daf5$t)
  
  
  Time <- c(0,meants1,meants2,meants3,meants4,meants5)
  std_dev <- c(0,sdts1,sdts2,sdts3,sdts4, sdts5)
  algo= rep(i,6)
  L<-c(0,200,500,1000,1500,2000)
  
  
  df4 <- data.frame(Time, std_dev, algo,L)
  df4$algo<-as.factor(df4$algo)
  
  df<-rbind(df,df4)
}
  

######################################################


colnames(df) <- c('Time', 'std_dev', 'Algorithm', 'L')

p <- ggplot(df, aes(L, Time)) + 
  geom_line(aes(colour=Algorithm)) + 
  geom_ribbon(aes(ymin=Time-(1.96/22.3)*std_dev, ymax=Time+(1.96/22.3)*std_dev, fill=Algorithm), alpha=.2) +
  theme_classic() +labs(x = "Number of influece samples", y="Expected size of the spread") +
  ggtitle(expression(Seed~set~k==1))+theme(plot.title = element_text(hjust = 0.5))
p


saving_name=paste("epsilon_",as.character(epsilon),"_k_",as.character(k), ".pdf",sep = "")

ggsave(saving_name)


