############################################################################################################
#Project: STAT 628 Module1 Body Fat Caculator
#Authro: Shuyi Qu      squ26@wisc.edu
#        Guoli Liu     gliu66@wisc.edu
#        Linhai Zhang  lzhang472@wisc.edu
#Date:   2018/2/5
#Note:   This R code uses a lot of packages, 
#        please be patient to install them. Thank you.
############################################################################################################


############################################################################################################
###Include the packages
library(dplyr)
library(ggplot2)
library(GGally)
library(reshape)
library(reshape2)
library(ggcorrplot)
library(ggpubr)

###Include the data
bodyfat=as_tibble(read.csv("BodyFat.csv"))
############################################################################################################


############################################################################################################
###Part 1 Data wrangling

##1.1 Detect if there is age effect

#People's general healthy condition may change when aging,
#So we use plots to see if there is any evidence that age
#has effect on other variables.

#p1.1.1
attach(bodyfat)
bodyfat$AgeRange=rep(0,252)
bodyfat$AgeRange[which(AGE<=36)]="<36"
bodyfat$AgeRange[which(AGE>36)]="36~45"
bodyfat$AgeRange[which(AGE>45)]="45~54"
bodyfat$AgeRange[which(AGE>54)]=">54"

bf1=ggboxplot(bodyfat,"AgeRange","BODYFAT",col="AgeRange")
ds1=ggboxplot(bodyfat,"AgeRange","DENSITY",col="AgeRange")
ag1=ggboxplot(bodyfat,"AgeRange","AGE",col="AgeRange")
we1=ggboxplot(bodyfat,"AgeRange","WEIGHT",col="AgeRange")
hi1=ggboxplot(bodyfat,"AgeRange","HEIGHT",col="AgeRange")

ggarrange(bf1,ds1,ag1,we1,hi1,common.legend=T,ncol=2,nrow=3)
detach(bodyfat)

#The plot shows that other variables seem to be uniform according to age,
#So no need to deal with age.

##1.2 Detect the odd points
summary(bodyfat)
#By the summary, we detect 5 odd points:
#42th  : HEIGHT is 29.5
#39ht  : WEIGHT is 363.15
#79th  : AGE is 81
#128th : BODYFAT is 0
#216th : BODYFAT is 45.1 (nearly impossible for human)

#Then we use a quantile table to investigate them
#p1.2.1
library(gridExtra)
library(grid)
attach(bodyfat)
od1=bodyfat[which(BODYFAT==min(BODYFAT)),]
od2=bodyfat[which(BODYFAT==max(BODYFAT)),]
od3=bodyfat[which(DENSITY==min(DENSITY)),]
od4=bodyfat[which(DENSITY==max(DENSITY)),]
od5=bodyfat[which(AGE==max(AGE)),]
od6=bodyfat[which(WEIGHT==max(WEIGHT)),]
od7=bodyfat[which(HEIGHT==min(HEIGHT)),]
odd=rbind.data.frame(od1,od2,od3,od4,od5,od6,od7)
odd=unique(odd)
odrank=matrix(0,ncol=7,nrow=5)
odrank[,1]=odd$IDNO
odrank[,7]=odd$AgeRange
for(i in 1:5){
  odrank[i,2]=round(ecdf(BODYFAT)(odd[i,2]),3)
  odrank[i,3]=round(ecdf(DENSITY)(odd[i,3]),3)
  odrank[i,4]=round(ecdf(AGE)(odd[i,4]),3)
  odrank[i,5]=round(ecdf(WEIGHT)(odd[i,5]),3)
  odrank[i,6]=round(ecdf(HEIGHT)(odd[i,6]),3)
}
odd1=as.data.frame(matrix(0,ncol=1,nrow=5))
odd1$IDNO=odrank[,1]
odd1$BODYFAT=odrank[,2]
odd1$DENSITY=odrank[,3]
odd1$AGE=odrank[,4]
odd1$WEIGHT=odrank[,5]
odd1$HEIGHT=odrank[,6]
odd1$AgeRange=odrank[,7]
odd1=odd1[,-1]
detach(bodyfat)
o1=tableGrob(odd1)
o2=tableGrob(odd[,c(1:6,18)])
grid.arrange(top=c("Percentile & Data Comparison"),o1,o2)

#Investigate if age has effect on correlations
#p1.2.2
bodyfat=as.data.frame(read.csv("BodyFat.csv"))
oddid=c(182,42)
bodyfat=bodyfat[-which(bodyfat$IDNO %in% oddid),]
attach(bodyfat)
bfdata=bodyfat[,-which(names(bodyfat) == "DENSITY")]
dsdata=bodyfat[,-which(names(bodyfat) == "BODYFAT")]

x=as.matrix(bodyfat[,-which(names(bodyfat) %in% c("BODYFAT","DENSITY"))])
y=as.matrix(bodyfat[,which(names(bodyfat) %in% c("BODYFAT","DENSITY"))])
corr=as.matrix(round(cor(x,y),1))

p.mat=as.matrix(cor_pmat(bodyfat))
p.mat=p.mat[-which(rownames(p.mat) %in% c("BODYFAT","DENSITY")),which(colnames(p.mat) %in% c("BODYFAT","DENSITY"))]

melted_corr=melt(corr)
ggheatmap=ggplot(melted_corr, aes(X1, X2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
ggheatmap+geom_text(aes(X1,X2,label=round(as.double(p.mat),2)),check_overlap = TRUE)+
  ggtitle("Correlation for full data")+
  labs(x="explanatory",y="BodyFat & Density")
detach(bodyfat)

#p1.2.3
bodyfat$AgeRange=as.logical(rep(0,250))
bodyfat$AgeRange[which(bodyfat$AGE<=45)]=FALSE
bodyfat$AgeRange[which(bodyfat$AGE>45)]=TRUE

attach(bodyfat)
x=as.matrix(bodyfat[which(AgeRange==TRUE),-which(names(bodyfat) %in% c("BODYFAT","DENSITY","AGE","AgeRange"))])
y=as.matrix(bodyfat[which(AgeRange==TRUE),which(names(bodyfat) %in% c("BODYFAT","DENSITY"))])
corr=as.matrix(round(cor(x,y),1))

p.mat=as.matrix(cor_pmat(bodyfat))
p.mat=p.mat[-which(rownames(p.mat) %in% c("BODYFAT","DENSITY","AGE","AgeRange")),which(colnames(p.mat) %in% c("BODYFAT","DENSITY"))]

melted_corr=melt(corr)
ggheatmap=ggplot(melted_corr, aes(X1, X2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
c1=ggheatmap+geom_text(aes(X1,X2,label=round(as.double(p.mat),2)),check_overlap = TRUE)+
  ggtitle("Correlation for Age under 45")+
  labs(x="explanatory",y="BodyFat & Density")

x=as.matrix(bodyfat[which(AgeRange !=TRUE),-which(names(bodyfat) %in% c("BODYFAT","DENSITY","AGE","AgeRange"))])
y=as.matrix(bodyfat[which(AgeRange !=TRUE),which(names(bodyfat) %in% c("BODYFAT","DENSITY"))])
corr=as.matrix(round(cor(x,y),1))

p.mat=as.matrix(cor_pmat(bodyfat))
p.mat=p.mat[-which(rownames(p.mat) %in% c("BODYFAT","DENSITY","AGE","AgeRange")),which(colnames(p.mat) %in% c("BODYFAT","DENSITY"))]

melted_corr=melt(corr)
ggheatmap=ggplot(melted_corr, aes(X1, X2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

c2=ggheatmap+geom_text(aes(X1,X2,label=round(as.double(p.mat),2)),check_overlap = TRUE)+
  ggtitle("Correlation for Age over 45")+
  labs(x="explanatory",y="BodyFat & Density")

ggarrange(c1,c2,nrow =2)

detach(bodyfat)

#From the quantile table we can see that
#the 42th, 128th and 216th data are odd,
#while 39th and 79th data are within the range.
#So we delete the 42th, 128th and 216th data.

bodyfat = bodyfat[-c(42,128,216),]
#change the unit of weight to kg, change the unit of height to meter
#save the cleaned data
bodyfat2 = bodyfat
bodyfat2$HEIGHT = bodyfat2$HEIGHT*0.0254
bodyfat2$WEIGHT = bodyfat2$WEIGHT*0.4536
bodyfat2$IDNO = 1:dim(bodyfat)[1]
write.csv(bodyfat2,file="BodyFat_cleaned.csv",row.names = F)
############################################################################################################


############################################################################################################
###Part 2 Model building

##2.1 
#Given the clean data set now, do variable selection to reduce the dimension of explanatory variables. 
#The Pearson Correlation and p-value give it an elementary sketch. 
d1=subset(bodyfat,bodyfat$AgeRange==T)
d2=subset(bodyfat,bodyfat$AgeRange==F)

s0=ggplot() + 
  geom_point(data = d1, aes(x = HIP, y = DENSITY), color = "red") +
  geom_point(data = d2, aes(x = HIP, y = DENSITY), color = "blue") +
  xlab('HIP') +
  ylab('DENSITY')+
  ggtitle("HIP VS DENSITY")
s1=ggplot() + 
  geom_point(data = d1, aes(x = ABDOMEN, y = DENSITY), color = "red") +
  geom_point(data = d2, aes(x = ABDOMEN, y = DENSITY), color = "blue") +
  xlab('ABDOMEN') +
  ylab('DENSITY')+
  ggtitle("ABDOMEN VS DENSITY")
s2=ggplot() + 
  geom_point(data = d1, aes(x = CHEST, y = DENSITY), color = "red") +
  geom_point(data = d2, aes(x = CHEST, y = DENSITY), color = "blue") +
  xlab('CHEST') +
  ylab('DENSITY')+
  ggtitle("CHEST VS DENSITY")

ggarrange(s0,s1,s2,nrow=3,ncol=1,common.legend = T)


##2.2
#This part contains two part:
#The frist part is to use stepwise regression with BIC to reduce the number of parameters.
#The second part is to try all combinations of the frist part's result to get the best model.

library(leaps)
data1=as.data.frame(bodyfat[,-which(names(bodyfat) == "DENSITY")])
data1=data1[-which(data1$IDNO %in% c(182,216,42)),]
attach(data1)
leaps=regsubsets(BODYFAT~.,data=data1,nbest=10)
# view results 
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.

library(MASS)
fit=lm(BODYFAT~.,data=data1)
step=stepAIC(fit, direction="both",k=log(250))
step$anova

library(DAAG)
library(partitions)
library(gtools)
library(car)
library(ModelMetrics)
xnames=c("WEIGHT","ABDOMEN","FOREARM","WRIST")
c_2=combinations(4,2)
output=list()
output$x=list()
output$mse=list()
output$vif=list()
output$bic=list()
output$inv=list()
output$adj.R2=list()
for(i in 1:dim(c_2)[1]){
  tempx=xnames[c_2[i,]]
  tempdata=data.frame(data1[,which(names(data1) %in% tempx)])
  tempdata$y=data1$BODYFAT
  attach(tempdata)
  mtemp=lm(y~.,data=tempdata)
  output$x[[i]]=tempx
  output$mse[[i]]=mse(tempdata$y,mtemp$fitted.values)
  output$vif[[i]]=round(vif(mtemp),2)
  output$bic[[i]]=BIC(mtemp)
  output$inv[[i]]=paste(round(quantile(mtemp$residuals,c(0.3,0.7)),1),collapse = " ~ ")
  output$adj.R2[[i]]=summary(mtemp)[[9]]
  detach(tempdata)
}

c_3=combinations(4,3)
output1=list()
output1$x=list()
output1$mse=list()
output1$vif=list()
output1$bic=list()
output1$inv=list()
output1$adj.R2=list()
for(i in 1:dim(c_3)[1]){
  tempx=xnames[c_3[i,]]
  tempdata=data.frame(bodyfat[,which(names(bodyfat) %in% tempx)])
  tempdata$y=bodyfat$BODYFAT
  attach(tempdata)
  mtemp=lm(y~.,data=tempdata)
  output1$x[[i]]=tempx
  output1$mse[[i]]=mse(tempdata$y,mtemp$fitted.values)
  output1$vif[[i]]=round(vif(mtemp),2)
  output1$bic[[i]]=BIC(mtemp)
  output1$inv[[i]]=paste(round(quantile(mtemp$residuals,c(0.3,0.7)),1),collapse = " ~ ")
  output1$adj.R2[[i]]=summary(mtemp)[[9]]
  detach(tempdata)
}

tempdata=data.frame(data1[,-which(names(data1) %in% c("CHEST","ANKLE","KNEE","DENSITY","BODYFAT"))])
tempdata$y=data1$BODYFAT
attach(tempdata)
mtemp=lm(y~.,data=tempdata)
stemp=step(mtemp,k=log(249))
typeof(step$anova)

library(gridExtra)
library(grid)
library(gtable)
taov=data.frame(matrix(0,ncol=6,nrow=13))
taov$STEP=step$anova[[1]]
taov$DF=step$anova[[2]]
taov$DEVIANCE=step$anova[[3]]
taov$RESIDUAL.DF=step$anova[[4]]
taov$RESIDUAL.DEVIANCE=step$anova[[5]]
taov$BIC=step$anova[[6]]
aov=tableGrob(taov[,7:12])
plot.new()

#The result of 1st part
grid.draw(aov)

x=c()
for(j in 1:dim(c_2)[1]){
  x[j]=paste(output$x[[j]],collapse="+")
}
for(k in 1:dim(c_3)[1]){
  x=c(x,paste(output1$x[[k]],collapse="+"))
}

mse=c()
for(j in 1:dim(c_2)[1]){
  mse[j]=output$mse[[j]]
}
for(k in 1:dim(c_3)[1]){
  mse=c(mse,output1$mse[[k]])
}
mse=round(mse,2)

vif=c()
for(j in 1:dim(c_2)[1]){
  vif[j]=paste(output$vif[[j]],collapse=" , ")
}
for(k in 1:dim(c_3)[1]){
  vif=c(vif,paste(output1$vif[[k]],collapse=" , "))
}

bic=c()
for(j in 1:dim(c_2)[1]){
  bic[j]=output$bic[[j]]
}
for(k in 1:dim(c_3)[1]){
  bic=c(bic,output$bic[[k]])
}
bic=round(bic,2)

inv=c()
for(j in 1:dim(c_2)[1]){
  inv[j]=output$inv[[j]]
}
for(k in 1:dim(c_3)[1]){
  inv=c(inv,output1$inv[[k]])
}

adj_r2=c()
for(j in 1:dim(c_2)[1]){
  adj_r2[j]=output$adj.R2[[j]]
}
for(k in 1:dim(c_3)[1]){
  adj_r2=c(adj_r2,output1$adj.R2[[k]])
}
adj_r2=round(adj_r2,2)
model_select=data.frame(cbind(x,mse,bic,vif,inv,adj_r2))
colnames(model_select)=c("Model","MSE","BIC","VIF","Resid 30~70% CI","Adj.R2")
model_select=model_select[order(model_select$BIC),]
model_select=tableGrob(model_select)
plot.new()

#The result of 2nd part
grid.draw(model_select)

#Overall, the model is bodyfit = abdomen + wrist
tempdata=data.frame(data1[,which(names(data1) %in% c("ABDOMEN","WRIST"))])
tempdata$y=data1$BODYFAT
attach(tempdata)
fit=lm(y~.,data=tempdata)
summary(fit)
############################################################################################################


############################################################################################################
###Part 3 Model diagnostics

##3.1 Outliers
library(car)
outlierTest(fit) # Bonferonni p-value for most extreme obs
qqPlot(fit, main="QQ Plot") #qq plot for studentized resid 

##3.2 leverage plots
leveragePlots(fit)

##3.3 Influential Observations
# Influential Observations
# added variable plots 
av.plots(fit)
# Cook's D plot
# identify D values > 4/(n-k-1) 
cutoff <- 4/((nrow(mtcars)-length(fit$coefficients)-2)) 
plot(fit, which=4, cook.levels=cutoff)
# Influence Plot 
influencePlot(fit,	id.method="identify", main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )

##3.4 Non-normality
# Normality of Residuals
# qq plot for studentized resid
qqPlot(fit, main="QQ Plot")
# distribution of studentized residuals
library(MASS)
sresid <- studres(fit) 
hist(sresid, freq=FALSE, 
     main="Distribution of Studentized Residuals")
xfit<-seq(min(sresid),max(sresid),length=40) 
yfit<-dnorm(xfit) 
lines(xfit, yfit)
############################################################################################################
