setwd("C:\\Kulas\\Research\\weighting")
library(ggplot2)  #data visualization

cond1_nnnn <- read.csv("norm_norm_norm_norm_summary.csv")
cond2_pppp <- read.csv("pos_pos_pos_pos_summary.csv")
cond3_gggg <- read.csv("neg_neg_neg_neg_summary.csv")
cond4_nnng <- read.csv("norm_norm_norm_neg_summary.csv")
cond5_nnpp <- read.csv("norm_norm_pos_pos_summary.csv")
cond6_pnng <- read.csv("pos_norm_norm_neg_summary.csv")
cond7_pgpg <- read.csv("pos_neg_pos_neg_summary.csv")
cond8_ppgg <- read.csv("pos_pos_neg_neg_summary.csv")

cond1_nnnn$cond <- "Cond1"
cond2_pppp$cond <- "Cond2"
cond3_gggg$cond <- "Cond3"
cond4_nnng$cond <- "Cond4"
cond5_nnpp$cond <- "Cond5"
cond6_pnng$cond <- "Cond6"
cond7_pgpg$cond <- "Cond7"
cond8_ppgg$cond <- "Cond8"


combo <- rbind(cond1_nnnn, cond2_pppp, cond3_gggg, cond4_nnng, cond5_nnpp, cond6_pnng, cond7_pgpg, cond8_ppgg)
names(combo)[names(combo)=="cond"] <- "Condition"


combo$cellrate.ma <- with(combo, SSmale*SSdepta)
combo$cellrate.mb <- with(combo, SSmale*SSdeptb)
combo$cellrate.fa <- with(combo, SSfemale*SSdepta)
combo$cellrate.gb <- with(combo, SSfemale*SSdeptb)
combo <- transform(combo, SD=apply(combo[18:21],1, sd, na.rm = TRUE))

combo$RR.all <- with(combo, PPmale*PPdepta*SSmale*SSdepta + PPmale*PPdeptb*SSmale*SSdeptb +
                       PPfemale*PPdepta*SSfemale*SSdepta + PPfemale*PPdeptb*SSfemale*SSdeptb)



combo$NS.ma <- with(combo, 10000*PPmale*PPdepta*SSmale*SSdepta)
combo$NS.mb <- with(combo, 10000*PPmale*PPdeptb*SSmale*SSdeptb)
combo$NS.fa <- with(combo, 10000*PPfemale*PPdepta*SSfemale*SSdepta)
combo$NS.fb <- with(combo, 10000*PPfemale*PPdeptb*SSfemale*SSdeptb)
combo$NS <- with(combo,NS.ma+NS.mb+NS.fa+NS.fb)

combo$cell.sma <- with(combo, NS.ma/NS)*100
combo$cell.smb <- with(combo, NS.mb/NS)*100
combo$cell.sfa <- with(combo, NS.fa/NS)*100
combo$cell.sfb <- with(combo, NS.fb/NS)*100
combo$cell.pma <- with(combo, PPmale*PPdepta)*100
combo$cell.pmb <- with(combo, PPmale*PPdeptb)*100
combo$cell.pfa <- with(combo, PPfemale*PPdepta)*100
combo$cell.pfb <- with(combo, PPfemale*PPdeptb)*100


#combo$PPma <- scale(combo$cell.pma)
#combo$PPmb <- scale(combo$cell.pmb)
#combo$PPfa <- scale(combo$cell.pfa)
#combo$PPfb <- scale(combo$cell.pfb)
#combo$TSma <- scale(combo$cell.sma)
#combo$TSmb <- scale(combo$cell.smb)
#combo$TSfa <- scale(combo$cell.sfa)
#combo$TSfb <- scale(combo$cell.sfb)
#combo$one     <- (combo$PPma - combo$TSma)^2
#combo$two     <- (combo$PPmb - combo$TSmb)^2
#combo$three   <- (combo$PPfa - combo$TSfa)^2
#combo$four    <- (combo$PPfb - combo$TSfb)^2

combo$one     <- (combo$cell.pma - combo$cell.sma)^2
combo$two     <- (combo$cell.pmb - combo$cell.smb)^2
combo$three   <- (combo$cell.pfa - combo$cell.sfa)^2
combo$four    <- (combo$cell.pfb - combo$cell.sfb)^2


#The value E is the chance expected median value of sum of squared differences and
#if the d is measured as differences in stand scores, we obtain Formular 2, E = 2k

combo$Chron.d2 <- with(combo, one+two+three+four)
combo$r.p      <- with(combo, (3.357*2-Chron.d2)/(3.357*2+Chron.d2))


summary(combo$r.p)
which(combo$r.p < -0.936, arr.ind=T) 


## RQ1(a):total response rate and misrepresentation.[RR]
## need to calculate and mark correlation on each case.
new.fig.1 <- ggplot(combo, aes(x=RR.all, y=PPSU)) + 
  geom_point(alpha=.3) + 
  xlim(0.36,0.81) +
  ylim(0, .3) + 
  xlab("Total Response Rate") + 
  ylab("Misrepresentation") +
  theme(panel.background = element_rect(fill = "transparent", colour = "#CCCCCC"),
        panel.grid.major.y = element_line(colour = "grey80")) + 
  scale_colour_gradient() +
  geom_smooth(method='loess', colour="red") +
  theme(axis.title.x = element_text(size = 8,colour = "black")) +
  theme(axis.text.x = element_text(size = 8, colour = "black")) +
  theme(axis.title.y = element_text(size = 8,colour = "black")) +  
  theme(axis.text.y = element_text(size = 8,colour = "black"))
new.fig.1 + facet_wrap(~Condition) +
  theme(strip.text.x = element_text(size = 10,colour= "black"))


## RQ1(b):nonresponse form and misrepresentation.[SD]
## need to calculate and mark correlation on each case.
new.fig.2 <- ggplot(combo, aes(x=SD, y=PPSU, color=SD)) + 
  geom_point(alpha=.15) + 
  xlim(0, .2) +
  ylim(0, .3) + 

  xlab("SD indicating Form of Nonresponse") + 
  ylab("Misrepresentation") +
  theme(panel.background = element_rect(fill = "transparent", colour = "#CCCCCC"),
        panel.grid.major.y = element_line(colour = "grey80")) + 
  scale_colour_gradient() +
  geom_smooth(method='loess', colour="red") +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8,colour = "black")) +
  theme(axis.text.x = element_text(size = 8, colour = "black")) +
  theme(axis.title.y = element_text(size = 8,colour = "black")) +
  theme(axis.text.y = element_text(size = 8,colour = "black"))
new.fig.2 + facet_wrap(~Condition) +
  theme(strip.text.x = element_text(size = 10,colour= "black"))


## extended RQ1(b): sample representativeness and misrepresentation.[rp]
## need to calculate and mark correlation on each case.
new.fig.3 <- ggplot(combo, aes(x=r.p, y=PPSU)) + 
  geom_point(alpha=.3) + 
  xlim(-1.0, 1.0) + 
  ylim(0, .3) + 
  xlab("Cattell's Profile Similarity Coefficient rp") + 
  ylab("Misrepresentation") +
  theme(panel.background = element_rect(fill = "transparent", colour = "#CCCCCC"),
        panel.grid.major.y = element_line(colour = "grey80")) + 
  scale_colour_gradient() +
  geom_smooth(method='loess', colour="red") +
  theme(axis.title.x = element_text(size = 8,colour = "black")) +
  theme(axis.text.x = element_text(size = 8, colour = "black")) +
  theme(axis.title.y = element_text(size = 8,colour = "black")) +  
  theme(axis.text.y = element_text(size = 8,colour = "black"))
new.fig.3 + facet_wrap(~Condition) +
  theme(strip.text.x = element_text(size = 10,colour= "black"))


## RQ2 what impact does the application of weights have on both biased and unbiased sample estimates
## original figure 1 to 3

