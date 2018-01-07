####ASSIGNMENT GROUP3####

##read in datasets
data <- read.delim("~/UNI/Master/Fächer Ausland/Behaviour Genetics/Assignment3/mouse_hsd_ko_processed_data.txt",header=TRUE)
sample_data <- read.delim("~/UNI/Master/Fächer Ausland/Behaviour Genetics/Assignment3/mouse_hsd_ko_sample_data.txt", header=TRUE)

##load library for inspection
library(ggplot2)
library(limma)

expr<-data[,2:9]
row.names(expr) <- data$X

##make principle components
pca <- prcomp(t(expr), scale = TRUE)

##extract the components
components <- data.frame(pca$x[, 1:2])

##plot the components (assign plot into value which can later be called)
pca_plot <- qplot(x = PC1, y = PC2, data = components)

##as plot:
qplot(x=PC1, y=PC2, data=components)

##new plot with points coloured by group membership
components$group <- sample_data$group
pca_plot<-qplot(x=PC1, y = PC2, colour=group, data = components)
pca_plot

##differential expression
##assign a modelmatrix
model_matrix <- model.matrix(~group,data=sample_data)

model<-lmFit(log2(expr),model_matrix)
model<-eBayes(model)

##extract top probesets that are differently expressed
##look at them at 10% FalseDiscoveryRate
DGenes <- topTable(model,number=Inf,p.value=0.1)
DGenes

##statistical graphics
##plot expression of control and KO for one gene
plot_data <- data.frame(expression=as.vector(t(expr["1449038_at",])), group=sample_data$group)
scatterplot <- qplot(x=group, y=expression, data=plot_data, main="11b-HSD1", geom="jitter")
scatterplot

##plot with every gene
volcano_data <- topTable(model, number= Inf)
volcano_plot <- qplot(x=logFC, y=-log10(adj.P.Val), data=volcano_data)
volcano_plot

##for colours
volcano_data$threshold <- as.factor(abs(volcano_data$adj.P.Val) < 0.1)
volcano_plot <- qplot(x=logFC, y=-log10(adj.P.Val), colour=threshold, data=volcano_data)
volcano_plot

##to write the genes into a file
#write.csv(DGenes, file="~/UNI/Master/Fächer Ausland/Behaviour Genetics/Assignment3/Genetable")
