
library(dplyr)
library(ggplot2)
library(tidyr)
library(ROCR)
library(glmnet)

setwd("/Users/sichenghao/Desktop/DCM-Processing/data/")
features<-read.csv("features5.csv")
data<-read.csv("data5.csv")

#Data Clearning
features$files<-as.character(features$files)
features<-features%>%separate(files,c("N1","patient_number","N2","annotation_id","N3"))
features$patient_number<-as.numeric(features$patient_number)
features$annotation_id<-as.numeric(features$annotation_id)
features<-features[,-c(1,2,4,6)]

df<-left_join(features,data,by = c("patient_number","annotation_id"))


res.method<-rep(0,4)
res.treatment.diff<-rep(0,4)
#Model

#x<-as.matrix(df[,3:12])
x<-as.matrix(df[,c(4,6,7)])
y<-as.double(df[,17])
model1<-glmnet(x, y, family = "binomial",alpha = 1, lambda = NULL)
p1 <- predict(model1, newx=x, type="response")
coef(model1)
plot(model1)
coef(model1, s=min(model1$lambda))

#2 predict treatment by only using image features
model2 <- glm(treatment ~ features.no..1+features.no..3+features.no..4 + 
                health, data = df,family = "binomial")
#By adding health as a predictive factor, Residual Deviance is lower
predict <- predict(model2, type = 'response')
table(df$treatment, predict > 0.5)
ROCRpred <- prediction(predict, df$treatment)
ROCRperf <- performance(ROCRpred, 'tpr','fpr')
plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2,1.7))

#3 predict outcome from confounder
model3<-lm(data = df,survival.index ~features.no..1+features.no..3+features.no..4+health)
summary(model3)
p <- predict(model3, newx=df[,c(4,6,7,15)], type="response")
#4 predict outcome by only using treatment(biased)
model4<-lm(data = df,survival.index~treatment)
summary(model4)

res.method[1]<-"NA"
res.treatment.diff[1]<-model4$coefficients[2]

#5 add image feature to adjust
model5<-lm(data = df,survival.index~treatment+features.no..1+features.no..3+features.no..4)
summary(model5)
res.method[2]<-"Reg"
res.treatment.diff[2]<-model5$coefficients[2]




model6<-lm(data = df,survival.index~treatment+malignancy+health)
summary(model6)

model7<-lm(data = df,malignancy~features.no..1+features.no..3+features.no..4)
summary(model7)


dataset<-df[,c(4,6,7)]
dataset$Y = df$survival.index
dataset$W = as.numeric(df$treatment)-1

condmean <- glm(formula= Y ~ ., 
                data=dataset)
tauhat1x <- dataset %>%
  mutate(W = 1) %>%
  predict(condmean, type="response", newdata=.) %>%
  as.numeric()
tauhat0x <- dataset %>%
  mutate(W = 0) %>%
  predict(condmean, type="response", newdata=.) %>%
  as.numeric()

# Propensity score 
p <-randomForest(formula= I(factor(W)) ~ . -Y, 
                 data=dataset,
                 ntree=100,
                 type="classification",
                 seed=12345) %>%
  predict(., type="prob") %>% .[,2] %>% as.numeric()

# Double robust estimator
w <- dataset %>% pull(W)
y <- dataset %>% pull(Y)
idx<-c(which(p==0),which(p==1))
p<-p[-idx]
w<-w[-idx]
y<-y[-idx]
tauhat1x<-tauhat1x[-idx]
tauhat0x<-tauhat0x[-idx]
est1 <- w*(y - tauhat1x)/p + (1-w)*(y - tauhat0x)/(1-p)
est2 <- tauhat1x - tauhat0x
tauhat_dr_rf <- mean(est1, na.rm = TRUE) + mean(est2)
tauhat_dr_rf

condmean <- glm(formula= Y ~ ., 
                data=dataset)
tauhat1x <- dataset %>%
  mutate(W = 1) %>%
  predict(condmean, type="response", newdata=.) %>%
  as.numeric()
tauhat0x <- dataset %>%
  mutate(W = 0) %>%
  predict(condmean, type="response", newdata=.) %>%
  as.numeric()
w1 <- dataset %>% pull(W)
y1 <- dataset %>% pull(Y)
est1 <- w1*(y1 - tauhat1x)/p1 + (1-w1)*(y1 - tauhat0x)/(1-p1)
est2 <- tauhat1x - tauhat0x
tauhat_dr <- mean(est1, na.rm = TRUE) + mean(est2)
tauhat_dr

res.method[3]<-"Truth"
res.treatment.diff[3]<-3

res.method[4]<-"DR"
res.treatment.diff[4]<-tauhat_dr


res.df<-data.frame(res.method,res.treatment.diff)

ggplot(data = res.df[c(1,2,4),])+
  geom_col(aes(x = res.method,y = res.treatment.diff,fill = res.method))+
  xlab("")+ylab("Treatment Effect Difference")




