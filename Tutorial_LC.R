###Script to determine lethal concentration

##Set your working directory							
#setwd("~/Documents/R-Tutorial_LC")

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
G1=DataAll[DataAll$Species=="Aaegypti" & 
             DataAll$Lineage=="Pop01" & 
             DataAll$Substance=="S1" ,]
G2=DataAll[DataAll$Species=="Aaegypti" & 
             DataAll$Lineage=="Pop02"  & 
             DataAll$Substance=="S1" ,]
G3=DataAll[DataAll$Species=="Aalbopictus" & 
             DataAll$Lineage=="Pop01"  & 
             DataAll$Substance=="S1"  ,]

## fit logistic regression model for multiple conditions						
G1.glm=glm(data=G1, 
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
rownames(G1.lc, do.NULL = FALSE)
rownames(G1.lc)<-c("G1-50", "G1-90", "G1-95")
G1.lc 

G2.ld <- dose.p(G2.glm, p=c(0.50, 0.90, 0.95))  # from MASS
G2.ci <- G2.ld + attr(G2.ld, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)
G2.lc <- exp(cbind(G2.ld, attr(G2.ld, "SE"), G2.ci[,1], G2.ci[,2]))
dimnames(G2.lc)[[2]] <- c("LC", "SE", "Lower","Upper")
rownames(G2.lc, do.NULL = FALSE)
rownames(G2.lc)<-c("G2-50", "G2-90", "G2-95")
G2.lc 

G3.ld <- dose.p(G3.glm, p=c(0.50, 0.90, 0.95))  # from MASS
G3.ci <- G3.ld + attr(G3.ld, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)
G3.lc <- exp(cbind(G3.ld, attr(G3.ld, "SE"), G3.ci[,1], G3.ci[,2]))
dimnames(G3.lc)[[2]] <- c("LC", "SE", "Lower","Upper")
rownames(G3.lc, do.NULL = FALSE)
rownames(G3.lc)<-c("G3-50", "G3-90", "G3-95")
G3.lc


LC_final<-rbind(G1.lc, G2.lc, G3.lc)
sink("NAME_LC.txt")
print(LC_final)
sink()


##Prepare graphic visualization - install multiple packages
library(readr)
library(tibble)
library(stats)
library(dplyr)
library(ggplot2)

## Create data to plot							
Conc=seq(0, 200, 0.1)							
G1.pred=predict(G1.glm, newdata=data.frame(Concentration=Conc), 
                type="response")
Conc=seq(0, 200, 0.1)							
G2.pred=predict(G2.glm, newdata=data.frame(Concentration=Conc), 
                type="response")	
Conc=seq(0, 200, 0.1)							
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
G1.sg <- ggplot(G1.ndata, aes(x = Concentration, y = Mortality)) + 
  scale_x_continuous(limits = c(0,200), expand = expansion(mult = c(0, .02))) +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0, .02))) +
  geom_line(col="steelblue", size=1) + 
  geom_point(data=G1, shape=16, aes(x = Concentration, y = Mortality), col="steelblue", size=2) +
  labs(x = 'Concentration (mg/L)', y = 'Mortality', title = "a. Aedes aegypti - Pop01") + 
  geom_segment(aes(x = 0, y = 0.50, xend = G1.lc [1], yend = 0.50), 
               col="steelblue", size=1,linetype = "dashed") +
  geom_segment(aes(x = G1.lc [1], y = 0, xend = G1.lc [1], yend = 0.50), 
               col="steelblue", size=1, linetype = "dashed") +
  geom_segment(aes(x = 0, y = 0.90, xend = G1.lc [2], yend = 0.90), 
               col="steelblue", size=1,linetype = "dashed") +
  geom_segment(aes(x = G1.lc [2], y = 0, xend = G1.lc [2], yend = 0.90), 
               col="steelblue", size=1, linetype = "dashed") +
  theme(plot.title = element_text(size = 10, color = "black")) +
  theme(axis.title.y = element_text(size = 8, angle = 90))+
  theme(axis.title.x = element_text(size = 8, angle = 00)) +
  theme_bw()
G1.sg <- G1.sg + geom_ribbon(data = G1.ndata, aes(ymin = right_lwr, ymax = right_upr),
                                   alpha = 0.1, fill="steelblue", lty = 3, col = "steelblue") 
G1.sg

G2.sg <- ggplot(G2.ndata, aes(x = Concentration, y = Mortality)) + 
  scale_x_continuous(limits = c(0,200), expand = expansion(mult = c(0, .02))) +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0, .02))) +
  geom_line(col="red", size=1) + 
  geom_point(data=G2, shape=17, aes(x = Concentration, y = Mortality), col="red", size=2) +
  labs(x = 'Concentration (mg/L)', y = 'Mortality', title = "b. Aedes aegypti - Pop02") + 
  geom_segment(aes(x = 0, y = 0.50, xend = G2.lc [1], yend = 0.50), 
               col="red", size=1,linetype = "dashed") +
  geom_segment(aes(x = G2.lc [1], y = 0, xend = G2.lc [1], yend = 0.50), 
               col="red", size=1, linetype = "dashed") +
  geom_segment(aes(x = 0, y = 0.90, xend = G2.lc [2], yend = 0.90), 
               col="red", size=1,linetype = "dashed") +
  geom_segment(aes(x = G2.lc [2], y = 0, xend = G2.lc [2], yend = 0.90), 
               col="red", size=1, linetype = "dashed") +
  theme(plot.title = element_text(size = 10, color = "black")) +
  theme(axis.title.y = element_text(size = 8, angle = 90))+
  theme(axis.title.x = element_text(size = 8, angle = 00)) +
  theme_bw()
G2.sg <- G2.sg + geom_ribbon(data = G2.ndata, aes(ymin = right_lwr, ymax = right_upr),
                             alpha = 0.1, fill="red", lty = 3, col = "red") 
G2.sg

G3.sg <- ggplot(G3.ndata, aes(x = Concentration, y = Mortality)) + 
  scale_x_continuous(limits = c(0,200), expand = expansion(mult = c(0, .02))) +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0, .02))) +
  geom_line(col="darkgray", size=1) + 
  geom_point(data=G3, shape=18, aes(x = Concentration, y = Mortality), col="darkgray", size=2) +
  labs(x = 'Concentration (mg/L)', y = 'Mortality', title = "c. Aedes albopictus - Pop01") + 
  geom_segment(aes(x = 0, y = 0.50, xend = G3.lc [1], yend = 0.50), 
               col="darkgray", size=1,linetype = "dashed") +
  geom_segment(aes(x = G3.lc [1], y = 0, xend = G3.lc [1], yend = 0.50), 
               col="darkgray", size=1, linetype = "dashed") +
  geom_segment(aes(x = 0, y = 0.90, xend = G3.lc [2], yend = 0.90), 
               col="darkgray", size=1,linetype = "dashed") +
  geom_segment(aes(x = G3.lc [2], y = 0, xend = G3.lc [2], yend = 0.90), 
               col="darkgray", size=1, linetype = "dashed") +
  theme(plot.title = element_text(size = 10, color = "black")) +
  theme(axis.title.y = element_text(size = 8, angle = 90))+
  theme(axis.title.x = element_text(size = 8, angle = 00)) +
  theme_bw()
G3.sg <- G3.sg + geom_ribbon(data = G3.ndata, aes(ymin = right_lwr, ymax = right_upr),
                             alpha = 0.1, fill="darkgray", lty = 3, col = "darkgray") 
G3.sg

##Plot type 2
G1.dflc<-as.data.frame(G1.lc)
G1.dflc$Group<-c("G1", "G1", "G1")
G1.dflc$Type<-c("LC50", "LC90", "LC95")
G2.dflc<-as.data.frame(G2.lc)
G2.dflc$Group<-c("G2", "G2", "G2")
G2.dflc$Type<-c("LC50", "LC90", "LC95")
G3.dflc<-as.data.frame(G3.lc)
G3.dflc$Group<-c("G3", "G3", "G3")
G3.dflc$Type<-c("LC50", "LC90", "LC95")
LC.dflc<-rbind(G1.dflc, G2.dflc, G3.dflc)
LC.graph<- ggplot(data=LC.dflc, aes(x=Type, y=LC, colour = Group, group = Group, shape = Group)) +
  geom_point(size = 3, position = position_dodge(width = 0.5))+
  scale_color_manual(values = c("steelblue", "red", "darkgray")) +
  scale_shape_manual(values = c(16,17,18)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.5), size=1, width=.4) + 
  coord_flip()  +  # flip coordinates (puts labels on y axis)
#  geom_hline(yintercept=100, lty=2)  + 
  xlim("LC50", "LC90", "LC95") + ylim(0, 50) +  # add a dotted line at x=1 after flip
  labs(x= " ",y = 'Estimate', title = "d. Lethal Concentration") +
  theme_bw()  # use a white background
LC.graph

##Export figure with multiple graphics
library(gridExtra)
Graph1 <- arrangeGrob(G1.sg, G2.sg, G3.sg, LC.graph, nrow=2) #generates g
ggsave("figure1.png", Graph1, units="mm", width=340, height=340, dpi=600)





