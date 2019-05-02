
library(dplyr)
library(ggplot2)

setwd("/Users/sichenghao/Desktop/DCM-Processing/data/")
annotation<-read.csv("annotation_df.csv")

data<-annotation%>%select(patient_number,annotation_id,malignancy)

#generate treatment
patient<-unique(data$patient_number)
#Randomly sample patient's health N(0,2)
#Make sure there are treatment for every malignancy
len<-length(patient)
set.seed(111)
health<-rnorm(len,0,1)
patient.df<-data.frame(patient_number = patient,health)
data<-left_join(data,patient.df)
data$treatment.score = data$health+data$malignancy
# ave.malignancy<-data%>%group_by(patient_number)%>%summarise(ave.malignancy = mean(malignancy))
# plot(density(ave.malignancy$ave.malignancy))
# patient.df<-left_join(patient.df,ave.malignancy)
# patient.df<-patient.df%>%mutate(treatment.score = health+ave.malignancy)
q50<-quantile(data$treatment.score,probs = 0.5)

data$treatment <- cut(data$treatment.score, breaks=c(-Inf, q50, Inf), labels=c("t1","t2"))
data%>%group_by(malignancy,treatment)%>%count()#check positivity


#Simulate treatment effect:
count<-data%>%group_by(treatment)%>%count()
t1.len<-count$n[1]
t2.len<-count$n[2]
set.seed(111)
t1.effect<-rnorm(t1.len,3,1)
set.seed(111)
t2.effect<-rnorm(t2.len,6,1)
treat.effect<-c(t1.effect,t2.effect)
data<-data[order(data$treatment),]
data$treat.effect<-treat.effect
data$survival.index<-data$treat.effect-3*data$malignancy+data$health+rnorm(length(data),0,1)

write.csv(data,file="data5.csv")

model<-lm(data = data,survival.index~treatment)
summary(model)

