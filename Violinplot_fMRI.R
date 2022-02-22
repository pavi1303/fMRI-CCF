library(ggplot2)
library(ggpubr)
library(ggsignif)
setwd("E:/My samsung/Work drive/Final_resampled_results")

data1 = readxl::read_xlsx("Reformatted_results_AUC_v2.xlsx",sheet = "violin_data")

ggplot(data1, aes(x=factor(Group), y=Correlation_coefficient, fill=Cognition_status)) + 
  geom_violin(trim = FALSE) +
  theme_bw() +
  geom_point(position=position_jitterdodge(0.03),dodge.width=0.5,aes(group=Cognition_status),alpha=0.2,show.legend = FALSE)+
  #geom_jitter(color="black",shape=21,aes(group=Disease_status),width=0.1)+
  geom_boxplot(width=0.06, color="black",position=position_dodge(0.9),alpha=0.5,lwd=1.0,outlier.shape=NA,show.legend = FALSE)+
  scale_fill_manual(values=c("deepskyblue1", "darkgoldenrod1"),name = "Cognition status", labels = c("Normal cognition","MCI"))+
  labs(y = "Correlation coefficient")+
  ggtitle("Violin plot - Group comparison")+
  scale_x_discrete(labels=c("A" = "Grp1", "B" = "Grp2"))+
  theme_bw()+
  theme(text=element_text(size=18,  family="serif"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 20, colour = "black"), axis.text.y = element_text(size = 18))+
  theme(legend.position = "bottom", legend.title = element_blank(),legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),legend.text = element_text(size = 16, colour = "black"))+
  ylim(-5,7)

