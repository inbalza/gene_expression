library("limma")
library("edgeR")
library("Glimma")
library("RColorBrewer") 
library("Rtsne")

# Read table from all_counts.txt
countData <- read.table(file = "all_counts.txt",
                        sep = "\t", 
                        header = TRUE, 
                        stringsAsFactors = FALSE)
# Check first 5 rows
countData[1:5,]

# If first column contains genes (not counts) --> use as rownames
rownames(countData) <- countData$X
# Check first 5 rows again
countData[1:5,]
# Remove filename column
countData <- countData[,-1]
# Check first 5 rows again
countData[1:5,]

# Read metadata
metaData <- read.table(file = "SraRunTable.csv",
                       sep = ",", 
                       header = TRUE, 
                       stringsAsFactors = TRUE)
# Show colnames
colnames(metaData)


# Create DGEList object (requires edgeR) which will contains count data
DGEobject <- DGEList(counts = countData)
# Check dimensions 
dim(DGEobject)

# Add SraRunTable info to DGEobject$samples
DGEobject$samples <- cbind(DGEobject$samples, metaData[,c("Run","Batch","data_analysis_id","source_name")])
DGEobject$samples


cols <- brewer.pal(4,"Set2")
sample.col <- c(rep(cols[1],5),rep(cols[2],3),rep(cols[3],5),rep(cols[4],3))

# Barplot library sizes
par(mar = c(8,6,4,4))
bp <- barplot(rev(DGEobject$samples$lib.size*1e-6),
              col = rev(sample.col),
              axisnames = FALSE,
              main = "Library sizes",
              xlab = "Library size (millions)",
              xlim = c(0,60),
              horiz = TRUE)
axis(2, labels = rev(rownames(DGEobject$samples)), 
     at = bp, las = 1, cex.axis = 0.8)
text(rev(round(DGEobject$samples$lib.size*1e-6,2)),bp, 
     labels = rev(round(DGEobject$samples$lib.size*1e-6,2)), pos=4)


################################################################################
### DATA PRE-PROCESSING: Transformations from the raw scale
################################################################################

cpm <- cpm(DGEobject)
log.cpm <- cpm(DGEobject, log = TRUE)

# Genes that are not expressed
no.samples <- nrow(DGEobject$samples)
table(rowSums(DGEobject$counts==0)==no.samples)

DGEobject$counts[rowSums(DGEobject$counts==0)==no.samples,]

################################################################################
### DATA PRE-PROCESSING: Reduce subset of genes
################################################################################

keep.exprs <- rowSums(cpm>1)>=3
DGEfiltered <- DGEobject[keep.exprs,, keep.lib.sizes = FALSE]
dim(DGEfiltered)
nsamples = ncol(DGEfiltered)

# Now transform subset to get filtered+transformed cpmF (and log.cpmF)
cpmF <- cpm(DGEfiltered)
log.cpmF <- cpm(DGEfiltered, log = TRUE)

# Density plots of log-CPM values for raw pre-filtered data
# and post-filtered data (yF1000R)
par(mfrow = c(1,2))

# Density plot part A: raw data pre-filtered
plot(density(log.cpm[,1]), 
     col = sample.col[1], lwd = 2, las = 1, ylim = c(0,0.25),
     main = "A. Raw data", xlab = "Log-cpm")
abline(v = 0, lty = 3)
for (i in 2:nsamples){
  den <- density(log.cpm[,i])
  lines(den$x, den$y, col = sample.col[i], lwd = 2)
}
legend("topright", rownames(DGEobject$samples), 
       text.col = sample.col, cex = 0.5, bty = "n")
# Density plot part B: post-filtered
plot(density(log.cpmF[,1]), 
     col = sample.col[1], lwd = 2, las = 2, ylim = c(0,0.25),
     main = "B. Filtered data", xlab = "Log-cpm")
abline(v = 0, lty = 3)
for (i in 2:nsamples){
  den <- density(log.cpmF[,i])
  lines(den$x, den$y, col = sample.col[i], lwd = 2)
}
legend("topright", rownames(DGEfiltered$samples), 
       text.col = sample.col, cex = 0.5, bty="n")

################################################################################
### DATA PRE-PROCESSING: Normalising gene expression distributions
################################################################################

DGEf.notnorm <- DGEfiltered 
DGEf.norm <- calcNormFactors(DGEfiltered, method = "TMM")
DGEf.norm$samples$norm.factors

par(mfrow = c(1,2))
par(mar = c(8,4,2,2))

# Boxplot unnormalised data
DGEf.notnorm$samples$norm.factors
log.cpmF.notnorm <- cpm(DGEf.notnorm, log = TRUE)
boxplot(log.cpmF.notnorm, las = 2, col = sample.col, cex = 0.9, 
        main = "A. Unnormalised data", ylab = "Log-cpm")

# Boxplot normalised data
DGEf.norm$samples$norm.factors
log.cpmF.norm <- cpm(DGEf.norm, log = TRUE)
boxplot(log.cpmF.norm, las = 2, col = sample.col, cex = 0.9, 
        main = "B. Normalised data", ylab = "Log-cpm")

# --> only when scaling factors are not close to 1 you will see a big difference

################################################################################
### DATA PRE-PROCESSING: Unsupervised clustering of samples
################################################################################

par(mfrow = c(1,2), mar = c(6,4,4,2),oma = c(0, 0, 2, 0))

plotMDS(log.cpmF.norm, 
        pch = 19,
        col = sample.col,
        cex = 1,
        # xlim = c(-4.0,4.0), 
        # ylim = c(-2.0,2.0),
        dim = c(1,2),
        xaxt = "n",
        yaxt = "n",
        main="A. Before correction")
axis(1,cex.axis=0.75)
axis(2,cex.axis=0.75)


# Interactive MDS plot using glMDSPlot from Glimma
glMDSPlot(log.cpmF.norm, 
          labels = DGEfiltered$samples$source_name,
          groups = metaData[,c("source_name","data_analysis_id","Run")], 
          launch = TRUE)

###############################################################
### Plot MDS with batch correction
###############################################################

library("sva")

my_batch <- DGEfiltered$samples$Batch
my_mod <- model.matrix(~as.factor(DGEfiltered$samples$source_name))

# use ComBat function to correct for different batches 
combat <- ComBat(dat=log.cpmF.norm, batch=my_batch, mod=my_mod)

plotMDS(combat, 
        pch = 19,
        col = sample.col,
        cex = 1, 
        dim = c(1,2),
        xaxt = "n",
        yaxt = "n",
        main="B. After correction")
axis(1,cex.axis=0.75)
axis(2,cex.axis=0.75)
mtext("MDS Plot before and after batch correction", outer = TRUE, cex = 1.5)

#glMDSPlot(combat, 
#          labels = DGEfiltered$samples$source_name,
#          groups = metaData[,c("source_name","data_analysis_id","Run")], 
#          launch = TRUE)



################################################################################
### DIFFERENTIAL EXPRESSION ANALYSIS: Design matrix and contrasts
################################################################################
colnames(DGEf.norm)

## DESIGN MATRIX
group <- c(rep("MO",5), rep("MY",3), rep("WO",5), rep("WY",3))

design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

## CONTRASTS for pairwise comparisons between groups 
## are set up in limma using the makeContrasts function.
contr.matrix <- makeContrasts(
  MOvsMY = MO-MY,
  MOvsWO = MO-WO,
  WOvsWY = WO-WY,
  MYvsWY = MY-WY,
  levels = colnames(design))
contr.matrix

################################################################################
### DIFFERENTIAL EXPRESSION ANALYSIS: Removing heteroscedasticity
################################################################################
## Removing heteroscedasticity from count data
## In limma, linear modelling is carried out 
## on the log-CPM values ... by the voom function.
par(mfrow = c(1,2), mar = c(8,6,4,2))
vDGEf <- voom(DGEf.norm, design, plot = TRUE)
vDGEf
## Fitting linear models for comparisons of interest
## Linear modelling in limma is carried out using 
## the lmFit and contrasts.fit functions 
vfitDGEf <- lmFit(vDGEf, design)
vfitDGEf <- contrasts.fit(vfitDGEf, contrasts = contr.matrix)
efitDGEf <- eBayes(vfitDGEf)

plotSA(efitDGEf, main = "Final model")

################################################################################
### DIFFERENTIAL EXPRESSION ANALYSIS: Examining the DE genes
################################################################################
### Examining the number of DE genes
## Significance is defined using an adjusted p-value cutoff that is set at 5% by default.
dtEBayes <- decideTests(efitDGEf)
summary(dtEBayes)


## Some studies require more than an adjusted p-value cut-off. 
## When testing requires genes to have a log-FC that is significantly greater than 1 
## (equivalent to a 2-fold difference between cell types on the original scale).
tfitDGEf <- treat(vfitDGEf, lfc = log2(1.5))
dtDGEf <- decideTests(tfitDGEf)
summary(dtDGEf)

par(mfrow = c(1,2), mar = c(2,2,2,2))
vennDiagram(dt[,1:3], circle.col=c("purple", "dark green", "pink"), cex = 0.8)
vennDiagram(dt[,2:4], circle.col=c("dark green","pink","blue"), cex = 0.8)


## topTable for results using eBayes 
# C1=MO.MY, C2=MO.WO
MO.MY <- topTable(efitDGEf, coef = 1, n = Inf)
MO.WO <- topTable(efitDGEf, coef = 2, n = Inf)

topgenesC1 <- head(MO.MY, n = 1000)
topgenesC2 <- head(MO.WO, n = 1000)

# Top genes up and down sorted by logFC
topgenesC1.ordered <- topgenesC1[order(topgenesC1$logFC),] 
topgenesC2.ordered <- topgenesC2[order(topgenesC2$logFC),] 

top10genesC1.ordered <- head(topgenesC1.ordered, n = 5)
top10genesC1.ordered <- rbind(tail(topgenesC1.ordered, n = 5),top10genesC1.ordered)
top10genesC1.ordered <- top10genesC1.ordered[order(top10genesC1.ordered$logFC),] 
top10genesC1.ordered

top10genesC2.ordered <- head(topgenesC2.ordered, n = 5)
top10genesC2.ordered <- rbind(tail(topgenesC2.ordered, n = 5),top10genesC2.ordered)
top10genesC2.ordered <- top10genesC2.ordered[order(top10genesC2.ordered$logFC),] 
top10genesC2.ordered

################################################################################
### DIFFERENTIAL EXPRESSION ANALYSIS: graphical representations of DE genes
################################################################################

# Volcano plots
par(mfrow = c(2,2), mar = c(4,4,2,2))

cont <- c("MO/MY","MO/WO","WO/WY","MY/WY")
for (i in 1:length(cont)){
  volcanoplot(tfitDGEf, coef = i, main = cont[i],
              xlim = c(-4,4), ylim = c(0,12),
              names = rownames(tfitDGEf$coefficients))
  abline(h = 1.3, lty = 2, col = "red", lwd = 1.5)
  abline(v = c(log2(1.5),-log2(1.5)), lty = 2, col = "red", lwd = 1.5)
  
}

## MD plot
plotMD(tfit, column = 1, status = NULL, 
       main = "Title", xlim = c(-2,15))


## Create your own color palette to use in heatmap
color.palette.3 <- colorRampPalette(c("#0d50b2", "white", "#c5081a"))(n = 200)

## Heatmap topTable genes

# Get names of genes to show
genes10topTableC1 <- rownames(top10genesC1.ordered) # top
genes10topTableC2 <- rownames(top10genesC2.ordered) # top

topgenesTableC1 <- rownames(topgenesC1.ordered)
topgenesTableC2 <- rownames(topgenesC2.ordered)

## topTreat for results using treat 
topTreat.C1 <- topTreat(tfitDGEf, coef = 1, n = Inf)
topTreatgenesC1 <- head(topTreat.C1, n = 10)
topTreatgenesC1.ordered <- topTreatgenesC1[order(topTreatgenesC1$logFC),] 
topTreatgenesC1.ordered

topTreat.C2 <- topTreat(tfitDGEf, coef = 2, n = Inf)
topTreatgenesC2 <- head(topTreat.C2, n = 10)
topTreatgenesC2.ordered <- topTreatgenesC2[order(topTreatgenesC2$logFC),] 
topTreatgenesC2.ordered


# Packages required to use heatmap and str_count
library("gplots")
library("stringr")

par(mfrow = c(1,1), mar = c(8,6,4,2))

# Create heatmap top genes for MOvsWO
heatmap.2(vDGEf$E[topgenesTableC1,],scale = "row",
          Rowv = FALSE, # do not reorder rows!
          dendrogram = "column",
          col = color.palette.3, 
          trace = "none", 
          margins = c(8,6),
          #cexRow = 0.9,
          cexCol = 0.9,
          #labRow = topgenesTableC1, 
          labCol = group,
          key = TRUE,
          density.info = "none",
          lhei = c(2,20),
          lwid = c(2,2)
          
)


# Create heatmap top genes for MOvsWO
heatmap.2(vDGEf$E[topgenesTableC2,],scale = "row",
          Rowv = FALSE, # do not reorder rows!
          dendrogram = "column",
          col = color.palette.3, 
          trace = "none", 
          margins = c(8,6),
          #cexRow = 0.9,
          cexCol = 0.9,
          #labRow = topgenesTableC2, 
          labCol = group,
          key = TRUE,
          density.info = "none",
          lhei = c(2,20),
          lwid = c(2,2)
          
)      


## Heatmap topTreat geneslibrary
genestopTreatC1 <- rownames(topTreatgenesC1.ordered)
heatmap.2(vDGEf$E[genestopTreatC1,], 
          scale = "row",
          labRow = genestopTreatC1, 
          labCol = str_c(colnames(DGEf.norm), group, sep = "-"),
          col = color.palette.3, 
          trace = "none", density.info = "none",
          Rowv = FALSE,
          margin = c(10,6), lmat = rbind(4:3,2:1), 
          lhei = c(2,6), lwid = c(1.5,4), 
          dendrogram = "column")

genestopTreatC2 <- rownames(topTreatgenesC2.ordered)
heatmap.2(vDGEf$E[genestopTreatC2,], 
          scale = "row",
          labRow = genestopTreatC2, 
          labCol = str_c(colnames(DGEf.norm), group, sep = "-"),
          col = color.palette.3, 
          trace = "none", density.info = "none",
          Rowv = FALSE, # do not reorder rows!
          margin = c(10,6), lmat = rbind(4:3,2:1), 
          lhei = c(2,6), lwid = c(1.5,4), 
          dendrogram = "column")


################################################################################
### Check some genes
################################################################################
## Check counts for some genes of interest
expTable <- vDGEf$E
par(mfrow = c(1,2), mar = c(4,6,2,2))

geneSymbol <- "Tlr2"
DGEobject$counts[geneSymbol,]
mp <- barplot(rev(expTable[geneSymbol,]), xlim = c(0,6), horiz = TRUE,
              col = rev(sample.col), axes = F, axisnames = F, main = geneSymbol)
axis(1)
axis(2, labels = rev(colnames(expTable)), at = mp, las=1, cex.axis=0.8)
title(main = geneSymbol)

geneSymbol <- "Chodl"
mp <- barplot(rev(expTable[geneSymbol,]), xlim = c(0,6), horiz = TRUE,
              col = rev(sample.col), axes = F, axisnames = F, main = geneSymbol)
axis(1)
axis(2, labels = rev(colnames(expTable)), at = mp, las=1, cex.axis=0.8)
title(main = geneSymbol)

