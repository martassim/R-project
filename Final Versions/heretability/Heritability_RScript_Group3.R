#####ASSIGNMENT 1 - Group3#####

##read in data
Data <- as.data.frame(read.table(file="~/UNI/Master/Fächer Ausland/Behaviour Genetics/Assignment1/group3_data.csv", sep=",", header=T))

##make sure the factors are seen as factors
Data$animal<-as.factor(Data$animal)
Data$Age.category<-as.factor(Data$Age.category)
Data$Sex<-as.factor(Data$Sex)
Data$PCA3<-as.numeric(Data$PCA3)

##inspection of Data
head(Data)

##read in ped
Ped<- as.data.frame(read.table("~/UNI/Master/Fächer Ausland/Behaviour Genetics/Assignment1/beagle_pedigree.csv", sep=",", header=T))
for (x in 1:3) Ped[, x] <- as.factor(Ped[, x])

##inspection of Ped
head(Ped)

##load package
library(MCMCglmm)

##define priors
prior1.1<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

##animal model
model1.1<-MCMCglmm(PCA3~1,random=~animal,pedigree=Ped,data=Data,prior=prior1.1)

##Sol: Posterior distribution of location effects (and cutpoints for ordinal models)
plot(model1.1$Sol)

##The posterior distributions of the variance components (VCV)
plot(model1.1$VCV)

##new model with better fits (longer burn-in etc.)
model1.1 <- MCMCglmm(PCA3 ~ 1, random = ~animal, pedigree = Ped,data=Data,nitt=65000,thin=50,burnin=15000,prior=prior1.1,verbose=FALSE)

##plot again
plot(model1.1$Sol)
plot(model1.1$VCV)

##calculate autocorrelation
##all lag values should ideally be near zero
autocorr(model1.1$VCV)

##obtain estimates of the additive genetic and residual variance by calculating the modes of the posterior distributions
posterior.mode(model1.1$VCV)
HPDinterval(model1.1$VCV)

##ESTIMATING HERITABILITY
##apply heritability formula to each sample of the posterior distribution
posterior.heritability1.1 <- model1.1$VCV[, "animal"]/(model1.1$VCV[,"animal"]+model1.1$VCV[,"units"])

posterior.mode(posterior.heritability1.1)
HPDinterval(posterior.heritability1.1, 0.95)

##plot posterior distribution of heritability estimate
plot(posterior.heritability1.1)

##ADDING FIXED EFFECTS
##modify fixed effect portion
model1.2 <- MCMCglmm(PCA3 ~ Sex, random = ~animal, pedigree = Ped, data=Data,prior=prior1.1,nitt=65000,thin=50,burnin=15000, verbose=FALSE)

##assess significance of sex as fixed effect
posterior.mode(model1.2$Sol[, "Sex2"])
HPDinterval(model1.2$Sol[, "Sex2"], 0.95)
posterior.mode(model1.2$VCV)

plot(model1.2$Sol)
plot(model1.2$VCV)

##SEX effect was previously in residual variances, but no effect -> no changes

posterior.heritability1.2 <- model1.2$VCV[, "animal"]/(model1.2$VCV[, "animal"]+model1.2$VCV[,"units"])

posterior.mode(posterior.heritability1.2)
HPDinterval(posterior.heritability1.2, 0.95)

plot(posterior.heritability1.2)

##ADDING RANDOM EFFECTS
##new prior  ##variance in PCA3 explained by animal and Age
prior1.3 <- list(G = list(G1 = list(V = 1, n = 0.002), G2 = list(V = 1,n=0.002)),R=list(V=1,n=0.002))

model1.3 <- MCMCglmm(PCA3 ~ Sex, random = ~animal + Age.category, pedigree = Ped,data=Data,nitt=65000,thin=50,burnin=15000,prior=prior1.3, verbose=FALSE)

posterior.mode(model1.3$VCV)

plot(model1.3$Sol)
plot(model1.3$VCV)

##Adding heritability to model1.3
posterior.heritability1.3 <- model1.3$VCV[, "animal"]/(model1.3$VCV[, "animal"]+model1.3$VCV[,"Age.category"]+model1.3$VCV[,"units"])

posterior.mode(posterior.heritability1.3)
HPDinterval(posterior.heritability1.3, 0.95)

plot(posterior.heritability1.3)

##Let's compare the estimates of heritability from each of models 1.2, 1.3:
posterior.mode(posterior.heritability1.1)
posterior.mode(posterior.heritability1.2)
posterior.mode(posterior.heritability1.3)

##TESTING SIGNIFICANCE OF VARIANCE COMPONENTS
model1.1$DIC
model1.2$DIC
model1.3$DIC

##Additional output (for mean values etc.)
summary(model1.1)
summary(model1.2)
summary(model1.3)