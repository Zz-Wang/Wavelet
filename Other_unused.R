#Bayesian homework
setwd("/Users/zhwaaa/Google Drive/R/Bayesian/")
#library(xlsx)
#data: Survival times of gastric cancer patients
#http://www.mayo.edu/research/documents/gastrichtml/doc-10027680
#variable1 > survival time
#variable2 > status: 1=death, 0=censored
#variable3 > treatment: 1=chemotherapy 2= combined chemotherapy/radiation

data <- read.csv("text", header = FALSE, sep = "")
colnames(data) <- c("time", "status", "treat")
boxplot(time~treat, data=data)
hist(data$time)

#basic function of Kaplan-Meier survival curves r
library(survival)
#creating survial object
data$SurvObj <- with(data, Surv(time, status == 1))
km <- survfit(SurvObj ~ treat, data = data, conf.type = "log-log")
plot(km)
plot(km, mark.time = T)

#p398
