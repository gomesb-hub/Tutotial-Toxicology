###Script to determine Resistance Ratio

##Set your working directory							
#setwd("~/Documents/R")

##Load your bioassays data into R - Define input name						
DataAll<-read.table("input_tutorial.txt", h=TRUE)							

## calculate mortality and add a new column to the data frame							
DataAll$Mortality=DataAll$Dead/(DataAll$Total)							

## Calculate 95% CIs around proportions following Newcombe and include them in original dataset					
CIs<-c(10)							
for(i in 1:length(DataAll$Dead)) { 				
  CIs<-append(CIs, prop.test(DataAll$Dead[i], DataAll$Total[i], correct=TRUE)$conf.int)							
}							
CIs<-CIs[-1]							
CIs <- matrix(CIs, ncol=2, byrow=TRUE)							
CIs							

colnames(CIs, do.NULL = FALSE)							
colnames(CIs) <- c("Lower","Upper")							

DataAll<-cbind(DataAll, CIs)							
DataAll						

## create multiple data sets						
G1Ref=DataAll[DataAll$Species=="Aaegypti" & 
             DataAll$Lineage=="Rockefeller" & 
             DataAll$Substance=="S1" ,]
G2=DataAll[DataAll$Species=="Aaegypti" & 
             DataAll$Lineage=="Pop01"  & 
             DataAll$Substance=="S1" ,]
G3=DataAll[DataAll$Species=="Aaegypti" & 
             DataAll$Lineage=="Pop02"  & 
             DataAll$Substance=="S1"  ,]

## fit logistic regression model for multiple conditions						
G1.glm=glm(data=G1Ref, 
            cbind(Dead,Alive)~log(Concentration), 
            family="binomial"(link="logit"))	
G2.glm=glm(data=G2, 
            cbind(Dead,Alive)~log(Concentration), 
            family="binomial"(link="logit"))
G3.glm=glm(data=G3, 
            cbind(Dead,Alive)~log(Concentration), 
            family="binomial"(link="logit"))

##Calculate LCs for 3 populations - Install MASS
library(MASS)
G1.ld <- dose.p(G1.glm, p=c(0.50, 0.90, 0.95))  # from MASS
G1.ci <- G1.ld + attr(G1.ld, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)
G1.lc <- exp(cbind(G1.ld, attr(G1.ld, "SE"), G1.ci[,1], G1.ci[,2]))
dimnames(G1.lc)[[2]] <- c("LC", "SE", "Lower","Upper")
G1.lc 

G2.ld <- dose.p(G2.glm, p=c(0.50, 0.90, 0.95))  # from MASS
G2.ci <- G2.ld + attr(G2.ld, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)
G2.lc <- exp(cbind(G2.ld, attr(G2.ld, "SE"), G2.ci[,1], G2.ci[,2]))
dimnames(G2.lc)[[2]] <- c("LC", "SE", "Lower","Upper")
G2.lc 

G3.ld <- dose.p(G3.glm, p=c(0.50, 0.90, 0.95))  # from MASS
G3.ci <- G3.ld + attr(G3.ld, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)
G3.lc <- exp(cbind(G3.ld, attr(G3.ld, "SE"), G3.ci[,1], G3.ci[,2]))
dimnames(G3.lc)[[2]] <- c("LC", "SE", "Lower","Upper")
G3.lc

##Calculate RRs for G2 and G3 with G1Ref as reference - Install pairwiseCI
library(pairwiseCI)
#RR for G2
G2.RR50<- MOVERR(theta0 = G1.lc[1], ci0 = c(G1.lc[7], G1.lc[10]), 
                 theta1 = G2.lc[1], ci1 = c(G2.lc[7], G2.lc[10]), 
                 alternative = "two.sided")
G2.RR50est<-do.call(rbind.data.frame, G2.RR50[2])
G2.RR50ci<-do.call(rbind.data.frame, G2.RR50[1])
G2.RR50data<-cbind(G2.RR50est[1],G2.RR50ci[1],G2.RR50ci[2])
dimnames(G2.RR50data)[[2]] <- c("Ratio", "Lower","Upper")


G2.RR90<- MOVERR(theta0 = G1.lc[2], ci0 = c(G1.lc[8], G1.lc[11]), 
                 theta1 = G2.lc[2], ci1 = c(G2.lc[8], G2.lc[11]), 
                 alternative = "two.sided")
G2.RR90est<-do.call(rbind.data.frame, G2.RR90[2])
G2.RR90ci<-do.call(rbind.data.frame, G2.RR90[1])
G2.RR90data<-cbind(G2.RR90est[1],G2.RR90ci[1],G2.RR90ci[2])
dimnames(G2.RR90data)[[2]] <- c("Ratio", "Lower","Upper")

G2.RR95<- MOVERR(theta0 = G1.lc[3], ci0 = c(G1.lc[9], G1.lc[12]), 
                 theta1 = G2.lc[3], ci1 = c(G2.lc[9], G2.lc[12]), 
                 alternative = "two.sided")
G2.RR95est<-do.call(rbind.data.frame, G2.RR95[2])
G2.RR95ci<-do.call(rbind.data.frame, G2.RR95[1])
G2.RR95data<-cbind(G2.RR95est[1],G2.RR95ci[1],G2.RR95ci[2])
dimnames(G2.RR95data)[[2]] <- c("Ratio","Lower","Upper")

G2.RRtable<-rbind(G2.RR50data, G2.RR90data, G2.RR95data)
rownames(G2.RRtable, do.NULL = FALSE)
rownames(G2.RRtable)<-c("G2-RR50", "G2-RR90", "G2-RR95")
G2.RRtable$Group<-c("G2", "G2", "G2")
G2.RRtable$Type<-c("RR50", "RR90", "RR95")

#RR for G3
G3.RR50<- MOVERR(theta0 = G1.lc[1], ci0 = c(G1.lc[7], G1.lc[10]), 
                 theta1 = G3.lc[1], ci1 = c(G3.lc[7], G3.lc[10]), 
                 alternative = "two.sided")
G3.RR50est<-do.call(rbind.data.frame, G3.RR50[2])
G3.RR50ci<-do.call(rbind.data.frame, G3.RR50[1])
G3.RR50data<-cbind(G3.RR50est[1],G3.RR50ci[1],G3.RR50ci[2])
dimnames(G3.RR50data)[[2]] <- c("Ratio", "Lower","Upper")

G3.RR90<- MOVERR(theta0 = G1.lc[2], ci0 = c(G1.lc[8], G1.lc[11]), 
                 theta1 = G3.lc[2], ci1 = c(G3.lc[8], G3.lc[11]), 
                 alternative = "two.sided")
G3.RR90est<-do.call(rbind.data.frame, G3.RR90[2])
G3.RR90ci<-do.call(rbind.data.frame, G3.RR90[1])
G3.RR90data<-cbind(G3.RR90est[1],G3.RR90ci[1],G3.RR90ci[2])
dimnames(G3.RR90data)[[2]] <- c("Ratio", "Lower","Upper")

G3.RR95<- MOVERR(theta0 = G1.lc[3], ci0 = c(G1.lc[9], G1.lc[12]), 
                 theta1 = G3.lc[3], ci1 = c(G3.lc[9], G3.lc[12]), 
                 alternative = "two.sided")
G3.RR95est<-do.call(rbind.data.frame, G3.RR95[2])
G3.RR95ci<-do.call(rbind.data.frame, G3.RR95[1])
G3.RR95data<-cbind(G3.RR95est[1],G3.RR95ci[1],G3.RR95ci[2])
dimnames(G3.RR95data)[[2]] <- c("Ratio", "Lower","Upper")

G3.RRtable<-rbind(G3.RR50data, G3.RR90data, G3.RR95data)
rownames(G3.RRtable, do.NULL = FALSE)
rownames(G3.RRtable)<-c("G3-RR50", "G3-RR90", "G3-RR95")
G3.RRtable$Group<-c("G3", "G3", "G3")
G3.RRtable$Type<-c("RR50", "RR90", "RR95")

RR_final<-rbind(G2.RRtable, G3.RRtable)
sink("NAME_RR.txt")
print(RR_final)
sink()

##Prepare graphic visualization - install multiple packages
library(readr)
library(tibble)
library(stats)
library(dplyr)
library(ggplot2)

## Create data to plot							
Conc=seq(1, 200, 0.1)							
G1.pred=predict(G1.glm, newdata=data.frame(Concentration=Conc), 
                type="response")
Conc=seq(1, 200, 0.1)							
G2.pred=predict(G2.glm, newdata=data.frame(Concentration=Conc), 
                type="response")	
Conc=seq(1, 200, 0.1)							
G3.pred=predict(G3.glm, newdata=data.frame(Concentration=Conc), 
                type="response")	


## some data to predict at: 100 values over the range of leafHeight
G1.ndata <- data.frame(Concentration=Conc, Mortality=G1.pred)
G2.ndata <- data.frame(Concentration=Conc, Mortality=G2.pred)
G3.ndata <- data.frame(Concentration=Conc, Mortality=G3.pred)
## grad the inverse link function
ilink1l <- family(G1.glm)$linkinv
## add fit and se.fit on the **link** scale
G1.ndata <- bind_cols(G1.ndata, 
                      setNames(as_tibble(predict(G1.glm, G1.ndata, se.fit = TRUE)[1:2]),
                               c('fit_link','se_link')))
## create the interval and backtransform
G1.ndata <- mutate(G1.ndata,
                   fit_resp  = ilink1l(fit_link),
                   right_upr = ilink1l(fit_link + (2 * se_link)),
                   right_lwr = ilink1l(fit_link - (2 * se_link)))

## grad the inverse link function
ilink2l <- family(G2.glm)$linkinv
## add fit and se.fit on the **link** scale
G2.ndata <- bind_cols(G2.ndata,
                      setNames(as_tibble(predict(G2.glm, G2.ndata, se.fit = TRUE)[1:2]),
                               c('fit_link','se_link')))
## create the interval and backtransform
G2.ndata <- mutate(G2.ndata,
                   fit_resp  = ilink2l(fit_link),
                   right_upr = ilink2l(fit_link + (2 * se_link)),
                   right_lwr = ilink2l(fit_link - (2 * se_link)))

## grad the inverse link function
ilink3l <- family(G3.glm)$linkinv
## add fit and se.fit on the **link** scale
G3.ndata <- bind_cols(G3.ndata, 
                      setNames(as_tibble(predict(G3.glm, G3.ndata, se.fit = TRUE)[1:2]),
                               c('fit_link','se_link')))
## create the interval and backtransform
G3.ndata <- mutate(G3.ndata,
                   fit_resp  = ilink3l(fit_link),
                   right_upr = ilink3l(fit_link + (2 * se_link)),
                   right_lwr = ilink3l(fit_link - (2 * se_link)))

##Plot type 1
G2.graph <- ggplot(G1.ndata, aes(x = Concentration, y = Mortality)) + 
  scale_x_continuous(trans='log10', expand = expansion(mult = c(0, .02))) +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0, .02))) +
  geom_line(col="steelblue", size=1) + 
  geom_point(data=G1Ref, shape=16, aes(x = Concentration, y = Mortality), col="steelblue", size=2) +
  geom_line(data=G2.ndata, aes(x = Concentration, y = Mortality), col="red", size=1) +
  geom_point(data=G2, shape=17, aes(x = Concentration, y = Mortality), col="red", size=2) +
  labs(x = 'Concentration (mg/L)', y = 'Mortality', title = "a. Aedes aegypti - Pop1") + 
  theme(plot.title = element_text(size = 10, color = "black")) +
  theme(axis.title.y = element_text(size = 8, angle = 90))+
  theme(axis.title.x = element_text(size = 8, angle = 00)) +
  theme_bw()
G2.graph <- G2.graph + geom_ribbon(data = G1.ndata, aes(ymin = right_lwr, ymax = right_upr),
                                   alpha = 0.1, fill="steelblue", lty = 3, col = "steelblue") +
  geom_ribbon(data = G2.ndata, aes(ymin = right_lwr, ymax = right_upr),
              alpha = 0.1, fill="red", lty = 3, col="red")
G2.graph

G3.graph <- ggplot(G1.ndata, aes(x = Concentration, y = Mortality)) + 
  scale_x_continuous(trans='log10', expand = expansion(mult = c(0, .02))) +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0, .02))) +
  geom_line(col="steelblue", size=1) + 
  geom_point(data=G1Ref, shape=16, aes(x = Concentration, y = Mortality), col="steelblue", size=2) +
  geom_line(data=G3.ndata, aes(x = Concentration, y = Mortality), col="darkred", size=1) +
  geom_point(data=G3, shape=18, aes(x = Concentration, y = Mortality), col="darkred", size=2) +
  labs(x = 'Concentration (mg/L)', y = 'Mortality', title = "b. Aedes aegypti - Pop2") + 
  theme(plot.title = element_text(size = 10, color = "black")) +
  theme(axis.title.y = element_text(size = 8, angle = 90))+
  theme(axis.title.x = element_text(size = 8, angle = 00)) +
  theme_bw()
G3.graph <- G3.graph + geom_ribbon(data = G1.ndata, aes(ymin = right_lwr, ymax = right_upr),
                                   alpha = 0.1, fill="steelblue", lty = 3, col = "steelblue") +
  geom_ribbon(data = G3.ndata, aes(ymin = right_lwr, ymax = right_upr),
              alpha = 0.1, fill="darkred", lty = 3, col="darkred")
G3.graph

##Plot type 2
RR.graph<- ggplot(data=RR_final, aes(x=Type, y=Ratio, colour = Group, group = Group, shape = Group)) +
  geom_point(size = 4, position = position_dodge(width = 0.5))+
  scale_color_manual(values = c("red","darkred")) +
  scale_shape_manual(values = c(17,18)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.5), size=.8, width=.2) + 
  coord_flip()  +  # flip coordinates (puts labels on y axis)
  geom_hline(yintercept=1, lty=2, lwd=1)  + 
  geom_hline(yintercept=5, lty=3, lwd=1)  + 
  xlim("RR50", "RR90", "RR95") + ylim(0, 6) +  # add a dotted line at x=1 after flip
  labs(x= " ",y = 'Estimate', title = "Resistance Ratio") +
  theme_bw()  # use a white background
RR.graph

##Export figure with multiple graphics
library(gridExtra)
Graph2 <- arrangeGrob(G2.graph, G3.graph, RR.graph, nrow=1) #generates 
ggsave("figure2.png", Graph2, units="mm", width=340, height=100, dpi=600)



