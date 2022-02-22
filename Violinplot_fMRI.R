library(ggplot2)
library(ggpubr)
library(ggsignif)
setwd("E:/LRCBH/Projects/COBRE/Results/Documents/Excel")

data1 = readxl::read_xlsx("Cobre_fluency_study_v3.xlsx",sheet = "violin_plot")

ggplot(data1, aes(x=factor(Group), y=Correlation_coefficient, fill=Cognition_status)) + 
  geom_violin(trim = FALSE) +
  theme_bw() +
  geom_point(position=position_jitterdodge(0.1),dodge.width=0.5,aes(group=Cognition_status),alpha=0.5,show.legend = FALSE)+
  #geom_jitter(color="black",shape=21,aes(group=Disease_status),width=0.1)+
  geom_boxplot(width=0.06, color="white",position=position_dodge(0.9),alpha=0.5,lwd=1.0,outlier.shape=NA,show.legend = FALSE)+
  scale_fill_manual(values=c("deepskyblue1", "darkgoldenrod1"),name = "Cognition status", labels = c("Normal cognition","MCI"))+
  labs(y = "Correlation coefficient")+
  ggtitle("Violin plot - Group comparison")+
  scale_x_discrete(labels=c("A" = "", "B" = ""))+
  theme_bw()+
  theme(text=element_text(size=18,  family="serif"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 20, colour = "black"), axis.text.y = element_text(size = 18))+
  theme(legend.position = "bottom", legend.title = element_blank(),legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),legend.text = element_text(size = 16, colour = "black"))+
  ylim(-0.6,0.6)

