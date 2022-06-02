##### LOAD PACKAGES #####
library(tidyverse)
library(nlme)
library(emmeans)
library(car)
library(MuMIn)
library(multcompView)
library(gdata)
##### DEFINE FUNCTIONS #####
#Connecting letter reports for models with more than one position
#niter is the number of groups
#leniter is the number of comparisons within each group
letterreport=function(emsummary,niter,leniter,sumname){
  for (i in 1:niter){
    refs=c(((i*leniter)-(leniter-1)):(i*leniter)) #Find where each group starts and ends
    P=emsummary$contrasts$p.value.B[refs]
    cont=gsub(" ","",emsummary$contrasts$contrast[refs])
    names(P)=cont
    writeLines(c(sumname,i,"\n"),sep=" ")
    working.comp=multcompLetters(P)
    print(working.comp[[1]])
    rm(refs,working.comp)
  }
}
# Means and multiple comparisons
multcomp.int=function(modname){
  # Get means
  mod.bp=emmeans(modname, pairwise~depth|pos,adjust="none")
  mod.tr=emmeans(modname, pairwise~pos|depth,adjust="none")
  mod.bp.s=summary(mod.bp)
  mod.tr.s=summary(mod.tr)
  # Find number of comparisons
  Ltot=length(mod.bp.s$contrasts$p.value)+length(mod.tr.s$contrasts$p.value)
  # Adjust p-values
  mod.bp.s$contrasts$p.value.B=p.adjust(mod.bp.s$contrasts$p.value,
                                        method='bonferroni',
                                        n = Ltot)
  mod.tr.s$contrasts$p.value.B=p.adjust(mod.tr.s$contrasts$p.value,
                                        method='bonferroni',
                                        n = Ltot)
  
  # Print things
  writeLines("\nEstimated Marginal Means & Pairwise Comparisons Within Body Position")
  print(mod.bp.s)
  writeLines("\nEstimated Marginal Means & Pairwise Comparisons Within Depth")
  print(mod.tr.s)
  list(mod.bp.s,mod.tr.s)
}

##### SPEED DATA - KINE SET #####
keep(letterreport,multcomp.int,sure=TRUE)
#Load data
working=read.csv("CleanDataset_KINE.csv") %>%
  mutate(across(where(is.character),as.factor),
         depth=fct_rev(as.factor(depth))) %>%
  select(fish,trial,depth,spd)

#Generate model and test assumptions (normality and equal variance of residuals)
mod=lme(spd~depth,random=~1|fish,data=working,method="REML")
res=resid(mod)
qplot(res,geom="histogram")
plot(mod$data$depth,res)
leveneTest(res,mod$data$depth)
summary(mod)
anova.lme(mod,type='marginal')
r.squaredGLMM(mod)

# Get estimated marginal means and pairwise comparisons
mod.sig.s=summary(emmeans(mod,pairwise~depth,adjust="none"))
mod.full=emmeans(mod,~depth)
Lsig=length(mod.sig.s$contrasts$p.value)
mod.sig.s$contrasts$p.value.B=p.adjust(mod.sig.s$contrasts$p.value,
                                       method='bonferroni',
                                       n = Lsig)
P=mod.sig.s$contrasts$p.value.B
cont=gsub(" ","",mod.sig.s$contrasts$contrast)
names(P)=cont
multcompLetters(P)

##### CURVATURE COEFFICIENT - KINE SET #####
keep(letterreport,multcomp.int,sure=TRUE)
#Load data
working=read.csv("CleanDataset_KINE.csv") %>%
  mutate(across(where(is.character),as.factor),
         depth=fct_rev(as.factor(depth))) %>%
  select(fish,trial,depth,curvcoef)

#Generate model and test assumptions (normality and equal variance of residuals)
mod=lme(curvcoef~depth,random=~1|fish,data=working,
        weights=varIdent(form=~1|depth),
        method="REML")
res=resid(mod)
qplot(res,geom="histogram")
plot(mod$data$depth,res)
leveneTest(res,mod$data$depth)
summary(mod)
anova.lme(mod,type='marginal')
r.squaredGLMM(mod)

# Get estimated marginal means and pairwise comparisons
mod.sig.s=summary(emmeans(mod,pairwise~depth,adjust="none"))
mod.full=emmeans(mod,~depth)
Lsig=length(mod.sig.s$contrasts$p.value)
mod.sig.s$contrasts$p.value.B=p.adjust(mod.sig.s$contrasts$p.value,
                                       method='bonferroni',
                                       n = Lsig)
P=mod.sig.s$contrasts$p.value.B
cont=gsub(" ","",mod.sig.s$contrasts$contrast)
names(P)=cont
multcompLetters(P)

##### SWING DISTANCE - KINE SET #####
keep(letterreport,multcomp.int,sure=TRUE)
#Load data
datain=read.csv("CleanDataset_KINE.csv") %>%
  mutate(across(where(is.character),as.factor),
         depth=fct_rev(as.factor(depth))) %>%
  select(fish,trial,depth,swd.H,swd.20,swd.40,swd.60,swd.80,swd.T,spd)
working=stack(datain[,1:9])
names(working)=c('swd','pos')
working$fish=rep(datain$fish,6)
working$trial=rep(datain$trial,6)
working$depth=rep(datain$depth,6)
working$spd=rep(datain$spd,6)
working=working[complete.cases(working), ]

#Generate model and test assumptions (normality and equal variance of residuals)
mod=lme(swd~spd+depth*pos,random=~1|fish,data=working,
        weights=varIdent(form=~1|pos),
        method="REML")
res=resid(mod)
qplot(res,geom="histogram")
plot(mod$data$depth,res)
plot(mod$data$pos,res)
leveneTest(res,mod$data$depth)
leveneTest(res,mod$data$pos)
summary(mod)
anova.lme(mod,type='marginal')
r.squaredGLMM(mod)

pairsum=multcomp.int(mod)
letterreport(pairsum[[1]],6,28,"Position") #niter=n(pos)=6; leniter=n(comps within pos)=28
letterreport(pairsum[[2]],8,15,"Treatment") #niter=n(depth)=8; leniter=c(comps within depth)=15

##### WAVE FREQUENCY - KINE SET #####
keep(letterreport,multcomp.int,sure=TRUE)
#Load data
datain=read.csv("CleanDataset_KINE.csv") %>%
  mutate(across(where(is.character),as.factor),
         depth=fct_rev(as.factor(depth))) %>%
  select(fish,trial,depth,frq.40,frq.60,frq.80,frq.T,spd)
working=stack(datain[,1:7])
names(working)=c('frq','pos')
working$fish=rep(datain$fish,4)
working$trial=rep(datain$trial,4)
working$depth=rep(datain$depth,4)
working$spd=rep(datain$spd,4)
working=working[complete.cases(working), ]

#Generate model and test assumptions (normality and equal variance of residuals)
mod=lme(frq~spd*depth+pos,random=~1|fish,data=working,
        weights=varIdent(form=~1|depth),
        method="REML")
res=resid(mod)
qplot(res,geom="histogram")
plot(mod$data$depth,res)
plot(mod$data$pos,res)
leveneTest(res,mod$data$depth)
leveneTest(res,mod$data$pos)
summary(mod)
anova.lme(mod,type='marginal')
r.squaredGLMM(mod)

pairsum=multcomp.int(mod)
letterreport(pairsum[[1]],4,28,"Position") #niter=n(pos)=4; leniter=n(comps within pos)=28
letterreport(pairsum[[2]],8,6,"Treatment") #niter=n(depth)=8; leniter=c(comps within depth)=6

##### FIN FREQUENCY - KINE SET #####
keep(letterreport,multcomp.int,sure=TRUE)
working=read.csv("CleanDataset_KINE.csv") %>%
  mutate(across(where(is.character),as.factor),
         depth=fct_rev(as.factor(depth))) %>%
  select(fish,trial,depth,frq.F,spd)

mod=lme(frq.F~spd*depth,random=~1|fish,data=working,
        weights=varIdent(form=~1|depth),
        method="REML")
res=resid(mod)
qplot(res,geom="histogram")
plot(mod$data$depth,res)
leveneTest(res,mod$data$depth)
summary(mod)
anova.lme(mod,type='marginal')
r.squaredGLMM(mod)

# Get estimated marginal means and pairwise comparisons
mod.sig.s=summary(emmeans(mod,pairwise~depth,adjust="none"))
mod.full=emmeans(mod,~depth)
Lsig=length(mod.sig.s$contrasts$p.value)
mod.sig.s$contrasts$p.value.B=p.adjust(mod.sig.s$contrasts$p.value,
                                       method='bonferroni',
                                       n = Lsig)
P=mod.sig.s$contrasts$p.value.B
cont=gsub(" ","",mod.sig.s$contrasts$contrast)
names(P)=cont
multcompLetters(P)
##### FIN ROM - KINE SET #####
keep(letterreport,multcomp.int,sure=TRUE)
#Load data
working=read.csv("CleanDataset_KINE.csv") %>%
  mutate(across(where(is.character),as.factor),
         depth=fct_rev(as.factor(depth))) %>%
  select(fish,trial,depth,rom)

#Generate model and test assumptions (normality and equal variance of residuals)
mod=lme(rom~depth,random=~1|fish,data=working,
        weights=varIdent(form=~1|depth),
        method="REML")
res=resid(mod)
qplot(res,geom="histogram")
plot(mod$data$depth,res)
leveneTest(res,mod$data$depth)
summary(mod)
anova.lme(mod,type='marginal')
r.squaredGLMM(mod)

# Get estimated marginal means and pairwise comparisons
mod.sig.s=summary(emmeans(mod,pairwise~depth,adjust="none"))
mod.full=emmeans(mod,~depth)
Lsig=length(mod.sig.s$contrasts$p.value)
mod.sig.s$contrasts$p.value.B=p.adjust(mod.sig.s$contrasts$p.value,
                                       method='bonferroni',
                                       n = Lsig)
P=mod.sig.s$contrasts$p.value.B
cont=gsub(" ","",mod.sig.s$contrasts$contrast)
names(P)=cont
multcompLetters(P)

##### NOSE ELEVATION - KINE SET #####
keep(letterreport,multcomp.int,sure=TRUE)
#Load data
working=read.csv("CleanDataset_KINE.csv") %>%
  mutate(across(where(is.character),as.factor),
         depth=fct_rev(as.factor(depth))) %>%
  select(fish,trial,depth,elav) %>%
  drop_na(.) %>%
  droplevels(.)

#Generate model and test assumptions (normality and equal variance of residuals)
mod=lme(elav~depth,random=~1|fish,data=working,
        method="REML")
res=resid(mod)
qplot(res,geom="histogram")
plot(mod$data$depth,res)
leveneTest(res,mod$data$depth)
summary(mod)
anova.lme(mod,type='marginal')
r.squaredGLMM(mod)

# Get estimated marginal means and pairwise comparisons
mod.sig.s=summary(emmeans(mod,pairwise~depth,adjust="none"))
mod.full=emmeans(mod,~depth)
Lsig=length(mod.sig.s$contrasts$p.value)
mod.sig.s$contrasts$p.value.B=p.adjust(mod.sig.s$contrasts$p.value,
                                       method='bonferroni',
                                       n = Lsig)
P=mod.sig.s$contrasts$p.value.B
cont=gsub(" ","",mod.sig.s$contrasts$contrast)
names(P)=cont
multcompLetters(P)

##### BODY EMG DUTY FACTOR - EMG SET #####
keep(letterreport,multcomp.int,sure=TRUE)
#Load data
datain=read.csv("CleanDataset_EMG.csv") %>%
  mutate(across(where(is.character),as.factor),
         depth=fct_rev(as.factor(depth))) %>%
  select(fish,trial,depth,emgdf1,emgdf2,emgdf3,emgdf4,emgdf5)
working=stack(datain)
names(working)=c('df','pos')
working$fish=rep(datain$fish,5)
working$trial=rep(datain$trial,5)
working$depth=rep(datain$depth,5)
working$spd=rep(datain$spd,5)
working=working[complete.cases(working), ]

#Generate model and test assumptions (normality and equal variance of residuals)
mod=lme(df~depth+pos,random=~1|fish,data=working,
        method="REML")
res=resid(mod)
qplot(res,geom="histogram")
plot(mod$data$depth,res)
plot(mod$data$pos,res)
leveneTest(res,mod$data$depth)
leveneTest(res,mod$data$pos)
summary(mod)
anova.lme(mod,type='marginal')
r.squaredGLMM(mod)

pairsum=multcomp.int(mod)
letterreport(pairsum[[1]],5,10,"Position")
letterreport(pairsum[[2]],5,10,"Treatment")

##### FIN EMG DUTY FACTOR - EMG SET #####
keep(letterreport,multcomp.int,sure=TRUE)
#Load data
working=read.csv("CleanDataset_EMG.csv") %>%
  mutate(across(where(is.character),as.factor),
         depth=fct_rev(as.factor(depth))) %>%
  select(fish,trial,depth,emgdfF) %>%
  drop_na(.)
#Generate model and test assumptions (normality and equal variance of residuals)
mod=lme(emgdfF~depth,random=~1|fish,data=working,
        method="REML")
res=resid(mod)
qplot(res,geom="histogram")
plot(mod$data$depth,res)
leveneTest(res,mod$data$depth)
summary(mod)
anova.lme(mod,type='marginal')
r.squaredGLMM(mod)

# Get estimated marginal means
mod.sig.s=summary(emmeans(mod,pairwise~depth,adjust="none"))

##### BODY RIA - EMG SET #####
keep(letterreport,multcomp.int,sure=TRUE)
#Load data
datain=read.csv("CleanDataset_EMG.csv") %>%
  mutate(across(where(is.character),as.factor),
         depth=fct_rev(as.factor(depth))) %>%
  select(fish,trial,depth,ria1,ria2,ria3,ria4,ria5)
working=stack(datain)
names(working)=c('ria','pos')
working$fish=rep(datain$fish,5)
working$trial=rep(datain$trial,5)
working$depth=rep(datain$depth,5)
working$spd=rep(datain$spd,5)
working=working[complete.cases(working), ]

#Generate model and test assumptions (normality and equal variance of residuals)
mod=lme(ria~depth*pos,random=~1|fish,data=working,
        weights=varIdent(form=~1|pos),
        method="REML")
res=resid(mod)
qplot(res,geom="histogram")
plot(mod$data$depth,res)
plot(mod$data$pos,res)
leveneTest(res,mod$data$depth)
leveneTest(res,mod$data$pos)
summary(mod)
anova.lme(mod,type='marginal')
r.squaredGLMM(mod)

pairsum=multcomp.int(mod)
letterreport(pairsum[[1]],5,10,"Position")
letterreport(pairsum[[2]],5,10,"Treatment")

##### FIN RIA - EMG SET #####
keep(letterreport,multcomp.int,sure=TRUE)
#Load data
working=read.csv("CleanDataset_EMG.csv") %>%
  mutate(across(where(is.character),as.factor),
         depth=fct_rev(as.factor(depth))) %>%
  select(fish,trial,depth,riaF,spd) %>%
  drop_na(.)
#Generate model and test assumptions (normality and equal variance of residuals)
mod=lme(riaF~spd+depth,random=~1|fish,data=working,
        method="REML")
res=resid(mod)
qplot(res,geom="histogram")
plot(mod$data$depth,res)
leveneTest(res,mod$data$depth)
summary(mod)
anova.lme(mod,type='marginal')
r.squaredGLMM(mod)

# Get estimated marginal means and pairwise comparisons
mod.sig.s=summary(emmeans(mod,pairwise~depth,adjust="none"))
mod.full=emmeans(mod,~depth)
Lsig=length(mod.sig.s$contrasts$p.value)
mod.sig.s$contrasts$p.value.B=p.adjust(mod.sig.s$contrasts$p.value,
                                       method='bonferroni',
                                       n = Lsig)
P=mod.sig.s$contrasts$p.value.B
cont=gsub(" ","",mod.sig.s$contrasts$contrast)
names(P)=cont
multcompLetters(P)
