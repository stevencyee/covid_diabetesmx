library(readr); library(tidyverse); library(survival); library(mediation); library(ggpubr); library(rms)
library(survminer); library(haven); library(rsq); library(ResourceSelection); library(ggsci);library(timereg)
## Database management ####
setwd("C:/Users/HP-PC/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/COVID-19, Diabetes and obesity")
setwd("/Users/nefoantonio/UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/OMAR YAXMEHEN BELLO CHAVOLLA - COVID-19, Diabetes and obesity")

covid <- read_csv("200419COVID19MEXICO.csv")

covid$id<-paste0(str_pad(covid$ENTIDAD_RES, 2,pad = "0"),str_pad(covid$MUNICIPIO_RES,3, pad="0"))
covid1<-covid[,c(14:15,18:30,32:35)]
covid<-covid[,-c(14:15,18:30,32:35)]
covid1[covid1==2]<-0
covid1[covid1==97]<-NA;covid1[covid1==98]<-NA;covid1[covid1==99]<-NA
covid<-as.data.frame(cbind(covid, covid1))
covid$EMBARAZO[is.na(covid$EMBARAZO)]<-0
covid$TIPO_PACIENTE[covid$TIPO_PACIENTE==1]<-0;covid$TIPO_PACIENTE[covid$TIPO_PACIENTE==2]<-1
covid$covid<-NULL; covid$covid[covid$RESULTADO==1]<-1;covid$covid[covid$RESULTADO!=1]<-0
covid$edad65<-NULL;covid$edad65[covid$EDAD>=65]<-1;covid$edad65[covid$EDAD<65]<-0
covid$edad40<-NULL;covid$edad40[covid$EDAD>=40]<-0;covid$edad40[covid$EDAD<40]<-1
covid$diabetes_40<-NULL;covid$diabetes_40[covid$DIABETES==1 & covid$edad40==1]<-1;covid$diabetes_40[covid$DIABETES!=1 & covid$edad40!=1]<-0
covid$Mortalidad<-NULL; covid$Mortalidad[is.na(covid$FECHA_DEF)]<-0;covid$Mortalidad[is.na(covid$FECHA_DEF)==FALSE]<-1
covid$FECHA_DEF[is.na(covid$FECHA_DEF)]<-as.Date(Sys.Date())
covid$FU_time<-as.numeric(as.Date(covid$FECHA_DEF)-as.Date(covid$FECHA_SINTOMAS))
covid$FU_time[covid$FU_time>30]<-30
covid$Latencia<-as.numeric(as.Date(covid$FECHA_INGRESO)-as.Date(covid$FECHA_SINTOMAS))
covid$comorb<-covid$DIABETES+covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_d<-covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_ob<-covid$DIABETES+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$diab_ob<-2*covid$DIABETES+covid$OBESIDAD


write_sav(covid, "covid_mx.sav")

#### Positivity model ####
covid0<-covid %>% filter(RESULTADO!=3)

m1<-glm(covid~DIABETES*edad40+edad65+SEXO+HIPERTENSION+OBESIDAD+OBESIDAD+RENAL_CRONICA+
          ASMA+EPOC+CARDIOVASCULAR+INMUSUPR+TABAQUISMO, data=covid0, family="binomial")
OR<-as.data.frame(cbind(exp(coef(m1)),exp(confint(m1))))
OR<-OR[-c(1),]
OR$Covariate<-c("Diabetes", "Age <40", "Age >65","Male sex","Hypertension", "Obesity", "CKD", "CVD",
                "Asthma", "COPD", "Immunosupression", "Smoking", "Diabetes*Age <40")
colnames(OR)<-c("OR", "ciLow", "ciHigh", "Covariate")
OR$breaks<-seq(1:nrow(OR))
p1 <- ggplot(OR, aes(x = OR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = OR$breaks, labels = OR$Covariate) +
  ylab("") +
  xlab("All tests, Odds ratio (OR, 95%CI)")


###Diabetes
covid_diab<- covid0 %>% filter(RESULTADO!=3,DIABETES==1)

m2<-glm(covid~OBESIDAD+edad40+edad65+SEXO+HIPERTENSION+RENAL_CRONICA+
          ASMA+EPOC+CARDIOVASCULAR+INMUSUPR+TABAQUISMO, data=covid_diab, family="binomial")

OR<-as.data.frame(cbind(exp(coef(m2)),exp(confint(m2))))
OR<-OR[-c(1),]
OR$Covariate<-c("Obesity", "Age <40","Age >65", "Male sex","Hypertension", "CKD", "CVD",
                "Asthma", "COPD", "Immunosupression", "Smoking")
colnames(OR)<-c("OR", "ciLow", "ciHigh", "Covariate")
OR$breaks<-seq(1:nrow(OR))
p2 <- ggplot(OR, aes(x = OR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = OR$breaks, labels = OR$Covariate) +
  ylab("") +
  xlab("Diabetes, Odds ratio (OR, 95%CI)")

#Obesity
covid_ob<- covid0 %>% filter(RESULTADO!=3,OBESIDAD==1)
m2<-glm(covid~DIABETES*edad40+edad65+SEXO+HIPERTENSION+RENAL_CRONICA+
          ASMA+EPOC+CARDIOVASCULAR+INMUSUPR+TABAQUISMO, data=covid_ob, family="binomial")
OR<-as.data.frame(cbind(exp(coef(m2)),exp(confint(m2))))
OR<-OR[-c(1),]
OR$Covariate<-c("Diabetes", "Age <40", "Age >65","Male sex","Hypertension", "CKD", "CVD",
                "Asthma", "COPD", "Immunosupression", "Smoking", "Diabetes*Age <40")
colnames(OR)<-c("OR", "ciLow", "ciHigh", "Covariate")
OR$breaks<-seq(1:nrow(OR))
p3 <- ggplot(OR, aes(x = OR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = OR$breaks, labels = OR$Covariate) +
  ylab("") +
  xlab("Obesity, Odds ratio (OR, 95%CI)")

f1<-ggarrange(p1,p2,p3, labels = c("A", "B", "C"), nrow=1, ncol=3)

ggsave(f1,filename = "SuppFigure1.png", 
       bg = "transparent",
       width = 33.3, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

#### Lethality model ####
#Cox models
covid1<-covid%>%filter(RESULTADO!=3)
m2<-coxph(Surv(FU_time, Mortalidad)~covid*OBESIDAD+edad65+edad40+SEXO+HIPERTENSION+DIABETES+RENAL_CRONICA+CARDIOVASCULAR+
            ASMA+EPOC+INMUSUPR+TABAQUISMO, data=covid1)
HR<-as.data.frame(cbind(exp(coef(m2)),exp(confint(m2))))
HR$Covariate<-c("COVID-19", "Obesidad", "Age >65", "Age <40","Male sex","Hypertension", "Diabetes", "CKD", "CVD",
                "Asthma", "COPD", "Immunosupression", "Smoking","COVID-19*Obesity")
colnames(HR)<-c("HR", "ciLow", "ciHigh", "Covariate")
HR$breaks<-seq(1:nrow(HR))
h1 <- ggplot(HR, aes(x = HR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = HR$breaks, labels = HR$Covariate) +
  ylab("") +
  xlab("COVID-19 vs. non-COVID-19 mortality, Hazard ratio (HR, 95%CI)")

##Positive cases
covid1<-covid%>%filter(RESULTADO==1)

##Interaction modeling
m2<-coxph(Surv(FU_time, Mortalidad)~DIABETES*edad40+OBESIDAD+SEXO, data=covid1)
HR<-as.data.frame(cbind(exp(coef(m2)),exp(confint(m2))))
HR
#Figure
m2<-coxph(Surv(FU_time, Mortalidad)~DIABETES*edad40+edad65+OBESIDAD+SEXO+HIPERTENSION+RENAL_CRONICA+CARDIOVASCULAR+
            ASMA+EPOC+INMUSUPR+TABAQUISMO, data=covid1)
HR<-as.data.frame(cbind(exp(coef(m2)),exp(confint(m2))))
HR$Covariate<-c("Diabetes", "Age <40","Age >65", "Obesity","Male sex","Hypertension", "CKD", "CVD",
                "Asthma", "COPD", "Immunosupression", "Smoking","Diabetes*Age <40")
colnames(HR)<-c("HR", "ciLow", "ciHigh", "Covariate")
HR$breaks<-seq(1:nrow(HR))
h2 <- ggplot(HR, aes(x = HR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = HR$breaks, labels = HR$Covariate) +
  ylab("") +
  xlab("COVID-19 lethality, Hazard ratio (HR, 95%CI)")

#Diabetes
covid1<-covid%>%filter(RESULTADO==1, DIABETES==1)
m3<-coxph(Surv(FU_time, Mortalidad)~OBESIDAD+edad40+edad65+SEXO+HIPERTENSION+RENAL_CRONICA+CARDIOVASCULAR+
            ASMA+EPOC+INMUSUPR+TABAQUISMO, data=covid1)
HR<-as.data.frame(cbind(exp(coef(m3)),exp(confint(m3))))
HR$Covariate<-c("Obesity", "Age >65","Age <40","Male sex","Hypertension", "CKD", "CVD",
                "Asthma", "COPD", "Immunosupression", "Smoking")
colnames(HR)<-c("HR", "ciLow", "ciHigh", "Covariate")
HR$breaks<-seq(1:nrow(HR))
h3 <- ggplot(HR, aes(x = HR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = HR$breaks, labels = HR$Covariate) +
  ylab("") +
  xlab("COVID-19 lethality in DM, Hazard ratio (HR, 95%CI)")

##Model for all COVID-19 vs. no COVID-19
covid1<-covid%>%filter(RESULTADO!=3, DIABETES==1)
m3<-coxph(Surv(FU_time, Mortalidad)~covid*OBESIDAD+edad40+SEXO+HIPERTENSION+RENAL_CRONICA+CARDIOVASCULAR+
            ASMA+EPOC+INMUSUPR+TABAQUISMO, data=covid1)
HR<-as.data.frame(cbind(exp(coef(m3)),exp(confint(m3))))
HR$Covariate<-c("COVID-19", "Obesity","Age <40","Male sex","Hypertension", "CKD", "CVD",
                "Asthma", "COPD", "Immunosupression", "Smoking", "COVID-19*Obesity")
colnames(HR)<-c("HR", "ciLow", "ciHigh", "Covariate")
HR$breaks<-seq(1:nrow(HR))
sh1 <- ggplot(HR, aes(x = HR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = HR$breaks, labels = HR$Covariate) +
  ylab("") +
  xlab("COVID-19 vs non-COVID-19 lethality in DM, Hazard ratio (HR, 95%CI)")


### Obesity
covid1<-covid%>%filter(RESULTADO==1, OBESIDAD==1)
m4<-coxph(Surv(FU_time, Mortalidad)~DIABETES*edad40+edad65+SEXO+HIPERTENSION+RENAL_CRONICA+CARDIOVASCULAR+
            ASMA+EPOC+INMUSUPR+TABAQUISMO, data=covid1)
HR<-as.data.frame(cbind(exp(coef(m4)),exp(confint(m4))))
HR$Covariate<-c("Diabetes", "Age <40","Age >65","Male sex","Hypertension", "CKD", "CVD",
                "Asthma", "COPD", "Immunosupression", "Smoking", "Diabetes*Age <40")
colnames(HR)<-c("HR", "ciLow", "ciHigh", "Covariate")
HR$breaks<-seq(1:nrow(HR))
h4 <- ggplot(HR, aes(x = HR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = HR$breaks, labels = HR$Covariate) +
  ylab("") +
  xlab("COVID-19 lethality in obesity, Hazard ratio (HR, 95%CI)")

f2<-ggarrange(h1,h2,h3,h4,labels = c("A", "B", "C", "D"), nrow=2, ncol=2)

ggsave(f2,filename = "Figure1.png", 
       bg = "transparent",
       width = 32.2, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

#COVID mortality compared to non-COVID

##Model for all COVID-19 vs. no COVID-19
covid1<-covid%>%filter(RESULTADO!=3, OBESIDAD==1)
m3<-coxph(Surv(FU_time, Mortalidad)~DIABETES+covid+edad40+edad65+SEXO+HIPERTENSION+RENAL_CRONICA+CARDIOVASCULAR+
            ASMA+EPOC+INMUSUPR+TABAQUISMO, data=covid1)
HR<-as.data.frame(cbind(exp(coef(m3)),exp(confint(m3))))
HR$Covariate<-c("Diabetes", "COVID-19","Age <40","Age >65","Male sex","Hypertension", "CKD", "CVD",
                "Asthma", "COPD", "Immunosupression", "Smoking")
colnames(HR)<-c("HR", "ciLow", "ciHigh", "Covariate")
HR$breaks<-seq(1:nrow(HR))
sh2 <- ggplot(HR, aes(x = HR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = HR$breaks, labels = HR$Covariate) +
  ylab("") +
  xlab("COVID-19 vs. non-COVID-19 lethality in obesity,Hazard ratio (HR, 95%CI)")

sf2<-ggarrange(sh1,sh2,labels = c("A", "B", "C", "D"), nrow=1, ncol=2)

ggsave(sf2,filename = "SuppFigure2.png", 
       bg = "transparent",
       width = 32.2, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

#### Other outcomes ####

## ICU
covid2<- covid%>%filter(RESULTADO==1)
m3<-glm(UCI~SEXO+DIABETES*edad40+OBESIDAD+edad65, data=covid2, family="binomial")
HR<-as.data.frame(cbind(exp(coef(m3)),exp(confint(m3))))
HR<-HR[-c(1),]
HR$Covariate<-c("Male sex", "Diabetes", "Age<40", "Obesity","Age >65",  "Diabetes*Age<40")
colnames(HR)<-c("HR", "ciLow", "ciHigh", "Covariate")
HR$breaks<-seq(1:nrow(HR))
o1 <- ggplot(HR, aes(x = HR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = HR$breaks, labels = HR$Covariate) +
  ylab("") +
  xlab("ICU admission, Odds ratio (OR, 95%CI)")

# Mechanical ventilation
m3<-glm(INTUBADO~SEXO+DIABETES*edad40+OBESIDAD+edad65, data=covid2, family="binomial")
HR<-as.data.frame(cbind(exp(coef(m3)),exp(confint(m3))))
HR<-HR[-c(1),]
HR$Covariate<-c("Male sex", "Diabetes", "Age <40", "Obesity","Age >65",  "Diabetes*Age<40")
colnames(HR)<-c("HR", "ciLow", "ciHigh", "Covariate")
HR$breaks<-seq(1:nrow(HR))
o2 <- ggplot(HR, aes(x = HR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = HR$breaks, labels = HR$Covariate) +
  ylab("") +
  xlab("Mechanical ventilation, Odds ratio (OR, 95%CI)")

# Hospitalization
m3<-glm(TIPO_PACIENTE~SEXO+DIABETES*edad40+OBESIDAD+edad65, data=covid2, family="binomial")
HR<-as.data.frame(cbind(exp(coef(m3)),exp(confint(m3))))
HR<-HR[-c(1),]
HR
HR$Covariate<-c("Male sex", "Diabetes", "Obesity", "Age <40","Age >65",  "Diabetes*Age<40")
colnames(HR)<-c("HR", "ciLow", "ciHigh", "Covariate")
HR$breaks<-seq(1:nrow(HR))
o3 <- ggplot(HR, aes(x = HR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = HR$breaks, labels = HR$Covariate) +
  ylab("") +
  xlab("Hospitalization risk, Odds ratio (OR, 95%CI)")

f4<-ggarrange(o1,o2,o3,labels = c("A", "B", "C"), nrow=1, ncol=3)

ggsave(f4,filename = "Figure3.png", 
       bg = "transparent",
       width = 35, 
       height = 15,
       units=c("cm"),
       dpi = 400,
       limitsize = FALSE)

#### Cox regression mediating effects ####

##All cases
covid1<-covid%>%filter(RESULTADO==1)%>%dplyr::select(Mortalidad, FU_time, DIABETES, OBESIDAD, EDAD, SEXO)
source("mediation_aalen_cox_ci_pval.R")
method="Cox"

set.seed(123)

glm1<-glm(DIABETES~OBESIDAD+EDAD+SEXO, data=covid1, family="binomial");summary(glm1)
m2<-coxph(Surv(FU_time, Mortalidad)~DIABETES+OBESIDAD+EDAD+SEXO, data=covid1);summary(m2)
lambdas<-m2$coefficients
lambdas.var<-m2$var

med<-mediation_ci1(lambdas[2], lambdas[1], lambdas.var[2,2], lambdas.var[2,1], lambdas.var[1,1],
              glm1$coef[2], vcov(glm1)[2,2], G=10^6, method=method)

med_prop<-1.756423/(1.756423+1.832215)*100; med_propci1<-1.491702/(1.491702+1.546204)*100;med_propci2<-2.102500 /(2.170272+2.102500 )*100

paste0("The effect of dibetes in mortality is ",round(med_prop,2),"%, 95%CI(",round(med_propci1,2),"-",round(med_propci2,2),") mediated by obesity.")


#### Kaplan Meier Analysis####

#Positive COVID-19#
covid_km<-covid%>%filter(RESULTADO==1)%>% 
  dplyr::select(Mortalidad,FU_time,DIABETES,OBESIDAD,HIPERTENSION,edad40,comorb_d)%>% drop_na()

comorb_d_REC<-NULL; comorb_d_REC[covid_km$comorb_d==0]<-0; comorb_d_REC[covid_km$comorb_d>=1 & covid_km$comorb_d<=3]<-1; 
comorb_d_REC[covid_km$comorb_d>=4]<-2; table(comorb_d_REC)
comorb_d_DIC<-NULL; comorb_d_DIC[covid_km$comorb_d==0]<-0; comorb_d_DIC[covid_km$comorb_d>=1]<-1

covid_km<-cbind(covid_km,comorb_d_REC,comorb_d_DIC)
covid_km$DIABETES<-as.factor(covid_km$DIABETES)
covid_km$OBESIDAD<-as.factor(covid_km$OBESIDAD)
covid_km$comorb_d_REC<-as.factor(covid_km$comorb_d_REC)
covid_km$comorb_d_DIC<-as.factor(covid_km$comorb_d_DIC)
covid_km$HIPERTENSION<-as.factor(covid_km$HIPERTENSION)

covid_km$DIABETES <- factor(covid_km$DIABETES,levels = c(0,1),labels = c("No-DM", "DM"))
covid_km$HIPERTENSION <- factor(covid_km$HIPERTENSION,levels = c(0,1),labels = c("No-HT", "HT"))
covid_km$comorb_d_REC<- factor(covid_km$comorb_d_REC,levels = c(0,1,2),labels = c("0", "1-3", "≥4"))
covid_km$comorb_d_DIC<- factor(covid_km$comorb_d_DIC,levels = c(0,1),labels = c("No-Comorb", "Comorbs"))

##Diabetes and comorbidities

cox_covid3<- survfit(Surv(FU_time, Mortalidad) ~ comorb_d_DIC + DIABETES, data = covid_km)
surv_pvalue(cox_covid3,data=covid_km,method = "n",test.for.trend =TRUE)
cox_covid_figure3<-ggsurvplot(cox_covid3, data = covid_km, size = 1,palette = "bw",conf.int = F,risk.table = T,pval = TRUE,ggtheme = theme_classic(),xlab="Time (Days)",
                              ylab="Survival probability",
                              legend.labs = c("No-DM & No-Comorb", 
                                              "DM & No-Comorb",
                                              "No-DM & Comorb",
                                              "DM & Comorb"),
                              xlim = c(0,30),ylim= c(0.8,1.0),break.y.by= c(0.05),break.x.by= c(5),pval.coord = c(0, 0.80))

cox_covid_figure3<-cox_covid_figure3 + theme_survminer(base_size = 10,
                                                       base_family = "Arial",
                                                       font.x = c(10, "plain" ), 
                                                       font.y = c(10, "plain"),
                                                       font.caption = c(10, "plain"), 
                                                       font.tickslab = c(10, "plain"))

#Obesity and diabetes

cox_covid4<- survfit(Surv(FU_time, Mortalidad) ~ OBESIDAD + DIABETES, data = covid_km)
surv_pvalue(cox_covid4,data=covid_km,method = "n",test.for.trend =TRUE)
cox_covid_figure4<-ggsurvplot(cox_covid4, data = covid_km,
                              size = 1,                 
                              palette = "bw",
                              conf.int = F,  
                              risk.table = T,
                              pval = TRUE,
                              pval.method = TRUE,
                              log.rank.weights="n",
                              ggtheme = theme_classic(),
                              xlab="Time (Days)",
                              ylab="Survival probability",
                              legend.labs = c("No-DM & No-Ob", 
                                              "DM & No-Ob",
                                              "No-DM & Ob",
                                              "DM & Ob"),
                              xlim = c(0,30),
                              ylim= c(0.8,1.0),
                              break.y.by= c(0.05),
                              break.x.by= c(5),
                              pval.coord = c(0, 0.80),
                              cumevents = F)

cox_covid_figure4 <-cox_covid_figure4 + theme_survminer(base_size = 10,
                                    base_family = "Arial",
                                    font.x = c(10, "plain" ), 
                                    font.y = c(10, "plain"),
                                    font.caption = c(10, "plain"), 
                                    font.tickslab = c(10, "plain"))

#Obesity and comorbidities

cox_covid5<- survfit(Surv(FU_time, Mortalidad) ~ OBESIDAD + comorb_d_DIC, data = covid_km)
surv_pvalue(cox_covid5,data=covid_km,method = "n",test.for.trend =TRUE)
cox_covid5
cox_covid_figure5<-ggsurvplot(cox_covid5, data = covid_km,
                              size = 1,                 
                              palette = "bw",
                              conf.int = F,  
                              risk.table = T,
                              pval = TRUE,
                              pval.method = TRUE,
                              log.rank.weights="n",
                              ggtheme = theme_classic(),
                              xlab="Time (Days)",
                              ylab="Survival probability",
                              legend.labs = c("No-Ob & No-Comorb", 
                                              "No-Ob & Comorb",
                                              "Ob & Comorb"),
                              xlim = c(0,30),
                              ylim= c(0.8,1.0),
                              break.y.by= c(0.05),
                              break.x.by= c(5),
                              pval.coord = c(0, 0.80),
                              cumevents = F)

cox_covid_figure5 <-cox_covid_figure5 + theme_survminer(base_size = 10,
                                                        base_family = "Arial",
                                                        font.x = c(10, "plain" ), 
                                                        font.y = c(10, "plain"),
                                                        font.caption = c(10, "plain"), 
                                                        font.tickslab = c(10, "plain"))

#Diabetes and Early onset
cox_covid6<- survfit(Surv(FU_time, Mortalidad) ~ DIABETES + edad40, data = covid_km)
surv_pvalue(cox_covid6,data=covid_km,method = "n",test.for.trend =TRUE)
cox_covid_figure6<-ggsurvplot(cox_covid6, data = covid_km,
                              size = 1,                 
                              palette = "bw",
                              conf.int = F,  
                              risk.table = T,
                              pval = TRUE,
                              pval.method = TRUE,
                              log.rank.weights="n",
                              ggtheme = theme_classic(),
                              xlab="Time (Days)",
                              ylab="Survival probability",
                              legend.labs = c("No-DM & Age <40", 
                                              "DM & Age <40",
                                              "No-DM & Age ≥40",
                                              "DM & Age ≥40"),
                              xlim = c(0,30),
                              ylim= c(0.8,1.0),
                              break.y.by= c(0.05),
                              break.x.by= c(5),
                              pval.coord = c(0, 0.80),
                              cumevents = F)

cox_covid_figure6 <-cox_covid_figure6 + theme_survminer(base_size = 10,
                                                        base_family = "Arial",
                                                        font.x = c(10, "plain" ), 
                                                        font.y = c(10, "plain"),
                                                        font.caption = c(10, "plain"), 
                                                        font.tickslab = c(10, "plain"))

cox_covid_figure6
figure4<-arrange_ggsurvplots(list(cox_covid_figure3,
          cox_covid_figure5,
          cox_covid_figure4,
          cox_covid_figure6),
          ncol=2,nrow =2,)

ggsave(file = "Figure3.png", 
       print(figure4),
       bg = "transparent",
       width = 50, 
       height = 27.7,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

#### Mechanistic COVID-19 lethality score ####

## GENERAL
covid1<-covid%>%filter(RESULTADO==1)
mod1<-coxph(Surv((FU_time), Mortalidad)~edad65+EMBARAZO+DIABETES*edad40+OBESIDAD+NEUMONIA+RENAL_CRONICA+
              EPOC+INMUSUPR, data=covid1)
summary(mod1)

points<-round(coef(mod1)/min(abs(coef(mod1))));points
covid1$score<-covid1$edad65*2+covid1$DIABETES+covid1$OBESIDAD*2+covid1$RENAL_CRONICA*3+covid1$EMBARAZO*6+
  covid1$NEUMONIA*7+covid1$EPOC*2+covid1$INMUSUPR*2+covid1$edad40*(-5)+covid1$DIABETES*covid1$edad40*4
table(covid1$score)

mod1_pts<-coxph(Surv(FU_time, Mortalidad)~score, data=covid1)
summary(mod1_pts)


covid1$score_cat<-NULL;covid1$score_cat[covid1$score<=0]<-0
covid1$score_cat[covid1$score>=1 & covid1$score<=3]<-1;covid1$score_cat[covid1$score>=4 & covid1$score<=6]<-2
covid1$score_cat[covid1$score>=7 & covid1$score<=9]<-3;covid1$score_cat[covid1$score>=10]<-4
covid1$score_cat<-as.factor(covid1$score_cat)
table(covid1$score_cat)
mod1_ptcat<-coxph(Surv(FU_time, Mortalidad)~score_cat, data=covid1)
summary(mod1_ptcat)

##KM of point score
cox_score<- survfit(Surv(FU_time, Mortalidad) ~ score_cat, data = covid1)
surv_pvalue(cox_score,data=covid1,method = "n",test.for.trend =TRUE)
cox_score_cat<-ggsurvplot(cox_score, data = covid1, size = 1,palette = "bw",conf.int = T,risk.table = T,pval = TRUE,ggtheme = theme_classic(),xlab="Time (Days)",
                                                             ylab="Survival probability",
                                                             legend.labs = c("Low-Risk", 
                                                                             "Mid-Risk",
                                                                             "Moderate-Risk",
                                                                             "High-Risk",
                                                                             "Very High Risk"),
                                                             xlim = c(0,30),
                                                             ylim= c(0.8,1.0),
                                                             break.y.by= c(0.05),
                                                             break.x.by= c(5),
                                                             pval.coord = c(0, 0.80))

cox_score1<-cox_score_cat + theme_survminer(base_size = 10,
                                                base_family = "Arial",
                                                font.x = c(10, "plain" ), 
                                                font.y = c(10, "plain"),
                                                font.caption = c(10, "plain"), 
                                                font.tickslab = c(10, "plain"))

ggsave(file = "Figure5.png", 
       print(cox_score1),
       bg = "transparent",
       width = 30, 
       height = 16.6,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

## Cross-validation
set.seed(123)
mod1_pts1<-cph(Surv(FU_time, Mortalidad)~score_cat, data=covid1, x=T, y=T)
mod1_pts1


validate(mod1_pts1, method="crossvalidation",B=10, bw=FALSE, rule="aic",
         type="residual", sls=.05, aics=0, force=NULL, estimates=TRUE,
         pr=FALSE, dxy=TRUE)

#### Figure with histograms ####
covid4<-covid1%>% dplyr::select(FECHA_SINTOMAS,diab_ob, Mortalidad, score)%>%drop_na()
covid4$diab_ob<-as.factor(covid4$diab_ob)
levels(covid4$diab_ob)<-c("NoDM-NoOb", "Ob-NoDM", "DM-NoOb", "DM-Ob")
covid4$Mortalidad<-as.factor(covid4$Mortalidad)
levels(covid4$Mortalidad)<-c("Non-lethal cases", "Lethal cases")
g1<-ggplot(covid4, aes(x=FECHA_SINTOMAS, fill=diab_ob))+
  geom_histogram(col="black", binwidth = 30)+
  ylab("New confirmed COVID-19 cases")+
  xlab("Symptom onset (date)")+
  theme_classic()+
  geom_vline(xintercept = 5, size = 1, colour = "#FF3721",linetype = "dashed")+
  facet_wrap(~Mortalidad, scales = "free")+
  labs(fill="Diabetes")

g2<-ggplot(covid4, aes(x=score, fill=Mortalidad))+
  geom_histogram(aes(y=..density..),col="black",binwidth = 1)+
  ylab("Density")+xlab("Mechanistic COVID-19 score")+
  theme_classic()+
  geom_vline(xintercept = 5, size = 1, colour = "#FF3721",linetype = "dashed")+
  labs(fill="Lethality")

fig4<-ggarrange(g1, g2, labels = c("A", "B"), ncol=1, nrow=2)


ggsave(fig4,filename = "Figure4.2.png", 
       bg = "transparent",
       width = 20, 
       height = 25,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)
