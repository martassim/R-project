#### ASSIGNMENT 2 ####

##load library
library(qtl)
#library(help = "graphics")
#data<-read.cross(format=c("csv"),dir="~/Dropbox/Applied Ethology and Animal Biology/3. Behaviour Genetics/Assignments/Assignment 2 - QTL Mapping",file="QTL_group_3.csv")
data<-read.cross(format=c("csv"),dir="~/UNI/Master/Fächer Ausland/Behaviour Genetics/Assignment2",file="QTL_group_3.csv")


summary(data)
plot(data)

##individual plot call
plotMissing(data)
plotMap(data)
plotPheno(data, pheno.col="id") ##inspect if normally distributed, if there are outlyer
plotPheno(data, pheno.col="w212", main = "Individuals' Weight distribution", xlab = "Weight (g)")
plotPheno(data, pheno.col="sex")
plotPheno(data, pheno.col="comb_g", main = "Comb weight Distribution", xlab = "Comb weight (g)")

#to plot the map with the markers' label names 
plotMap(data, chr=c(3,5),show.marker.names=TRUE)

##missing plot rearranged according to the phenotype's value of the individuals
plotMissing(data, reorder=TRUE)
plotMissing(data, chr = 3, reorder = TRUE) #just one chromosome

#to drop the markers without genotypic information (in this data we have information for all the markers)
data <- drop.nullmarkers(data)
totmar(data)

##have 2 plots next to each other to see the differences in one glance
par(mfrow=c(1,2)) #the numbers define the proportions of the graph. If smaller, two graphs might fit in one window
plotMissing(data)
plotMissing(data,reorder=TRUE)
par(mfrow=c(1,1)) #to make graphs bigger so that just one fits in a window

##est.fr() estimates the sex-averaged recombination fraction between all pairs of genetic markers
data <- est.rf(data)

##plot recombination fractions
##red (yellow) corresponds to a small recombination fraction or a big LOD score, while blue is the reverse.
##Gray indicates missing values
plotRF(data)
plotRF(data, chr=5) #just for one chromosome

##re-estimate the genetic map with another error probability
##measure recombination events and from this the distances between the markers
newmap <- est.map(data, error.prob=0.01)
plotMap(data, newmap, ch=3)
plotMap(newmap, chr=c(3,5),show.marker.names = FALSE)
plotMap(newmap, show.marker.names = TRUE)
plotMap(data, newmap, show.marker.names=TRUE)

##to replace the oldmap with the new map. Not always convenient
#fake.f2 <- replace.map(fake.f2, newmap)


##IDENTIFYING GENOTYPING ERRORS
##calculate error LOD scores (the log of the probability that 2 genes are linked)
##LOD scores are calculated for each individual at each marker
##large scores indicate likely genotyping errors
data <- calc.errorlod(data, error.prob=0.01)

##list of genotypes that may be in error -> top.errorlod(cross, chr, cutoff=4, msg=TRUE)
top.error <- top.errorlod(data)
summary(top.error)
table(top.error$marker)

##inspect observed genotypes for a chromosome with likely genotyping errors
##Note that white = AA and black = AB (for a backcross)
##ind=c(1:10) to specify the individuals that we want to see
plotGeno(data, chr=3)
plotGeno(data, chr=3, ind=c(1:10), ylab = "Individual", main = "Genotyping error in chromosome 3")#, ylab = "Location of Markers (cM)")
plotGeno(data, chr=5, ind=c(1:10))

##to find the label of a marker and to delete it (index = the ordinal position of the marker in the chromosome)
find.marker(data, chr = 3, index = 3)
data<-drop.markers(data,"Gg_rs15777012")

newmap2 <- est.map(data, error.prob=0.01)
plotMap(newmap2, show.marker.names = FALSE)
data <- replace.map(data, newmap2)

##measurement of the proportion of missing genotype information in the genotype data. Proportion of markers genotyped in the study population
plotInfo(data)
plotInfo(data, chr=c(3,5))
plotInfo(data, chr=c(3,5), method="entropy")
plotInfo(data, chr=c(3,5), method="variance")
plotInfo(data, chr=c(3,5), method="both")

##QTLmapping. The argument step indicates the step size (in cM) at which the probabilities are calculated, and determines the step size at which later LOD scores are calculated
data <-calc.genoprob(data, step=1, error.prob=0.01)

#to perform a single QTL genome scan
out.em <- scanone(data, pheno.col ="comb_g") #using the EM algorithm ???
out.hk <- scanone(data, pheno.col ="comb_g", method="hk") #using the HK regression

#t-test to check if the means of the weight are different between sexes
t.test(comb_g ~ sex, data=sex_w212)
cor.test... #acabar

##pull phenotypes from data and pull them into one combined vector
sex <- pull.pheno(data, "sex")
w212 <- pull.pheno(data, "w212")
sex_w212 <- cbind(sex, w212)

##calculates QTL genotype probabilities, conditional on the available marker data
##The argument step indicates the step size (in cM) at which the probabilities are calculated, and determines the step size at which later LOD scores are calculated.
data <-calc.genoprob(data, step=1, error.prob=0.01)

##single-QTL genome scan with a normal model
out.em <- scanone(data, addcovar = sex_w212, pheno.col = "comb_g") ##MaximumLikelihood
out.hk <- scanone(data, method="hk", addcovar=sex_w212, pheno.col ="comb_g") ##Haley-Knott regression

##multiple imputation method
data <- sim.geno(data, step=2, n.draws=16, error.prob=0.01)
out.imp <- scanone(data, addcovar=sex_w212, pheno.col = "comb_g", method="imp") #Genome scan with a single QTL model, with possible allowance for covariates, using any of several possible models for the phenotype and any of several possible numerical methods.

summary(out.em)
summary(out.em, threshold=3)
summary(out.hk, threshold=3)
summary(out.imp, threshold=3)

##plot the results (up to 3 genome scans)
##high lotscore: QTL most likely 
plot(out.em, chr=c(3,5)) ##plot of the scan with a single QTL model 
plot(out.em, out.hk, out.imp, chr=c(3,5)) ##plot all three models
plot(out.hk, chr=c(3,5), col="blue", add=TRUE) 
plot(out.imp, chr=c(3,5), col="red", add=TRUE)#add=TRUE adds this plot to the last one
abline(h=3.08, lty="dotted", lwd=1, col ="black") #to add a dotted line in the threshold
##permutation test to get a genome-wide LOD significance threshold
operm.hk <- scanone(data, addcovar=sex_w212,pheno.col ="comb_g", method="hk",n.perm=1000)
summary(operm.hk, alpha=0.05)

##add permutation results in summary.scanone
##estimated genome-scan- adjusted p-values for inferred QTL
##get a report of all chromosomes meeting a certain significance level, with the corresponding LOD threshold calculated automatically
summary(out.hk, perms=operm.hk, alpha=0.05, pvalues=TRUE)

##fit specific model
##pull out data on fixed QTL locations
qc <- c(3, 5) ##get chromosome 1 and 4
qp <- c(15.8, 361.0) ##get positions
qtl <- makeqtl(data, chr=qc, pos=qp)

##create a formula which indicates which QTL are to be included in the fit
myformula <- y~Q1+Q2

##fit the model
out.fq <- fitqtl(data, pheno.col ="comb_g", qtl=qtl, formula = myformula, covar=sex_w212)
summary(out.fq)

#To check whether the QTLs found are likely to be right 
##get estimated effects of QTL
##The estimated effects are the differences between the heterozygote and homozygote groups
##The interaction effect is the difference between the differences
out.fq <- fitqtl(data,pheno.col ="comb_g", qtl=qtl, formula = myformula, get.ests=TRUE)
summary(out.fq)

##refine the estimated positions of the QTL
##does QTL estimation at the same time together -> maybe slight changes because it fits somewhere else better together/ explains more
##A QTL object may be provided, or one may specify the chromosomes and positions
revqtl <- refineqtl(data, qtl=qtl,formula = myformula, covar=sex_w212, pheno.col ="comb_g") ##copied until here
revqtl ##short summary ##a bit of the QTL moved, but not much
plot(revqtl)

##re-run fitqtl to get a fit with the new positions
##overall LOD scores should have increased slightly
out.fq2 <- fitqtl(data, qtl=revqtl,formula=myformula, covar=sex_w212, pheno.col ="comb_g")
print(summary(out.fq2))

##to know which allele does what/ how the effects are
find.marker() ##type position of QTL to get the nearest makrer and its name
effectplot(data,mname1="Gg_rs15749695", pheno.col = "comb_g")
effectplot(data,mname1="c5.loc360", pheno.col = "comb_g")
effectplot(data,mname1="c5.loc360",mname2="Gg_rs15749695", pheno.col = "comb_g")
