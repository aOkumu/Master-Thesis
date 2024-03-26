setwd("G:/My Drive/UHASSELT/ACADEMIC WORK/Second year/semister 2/Master Thesis Bionformatics-ICP")

# Loading necessary Libraries
library(SingleCellExperiment)
library(Biostrings)
library(RColorBrewer) 
library(ggplot2) 
library(dplyr)
library(phylotools)
library(scater)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(limma)
library(muscat)
library(purrr)
library(scran)
library(gridExtra)
library(glmmTMB)

#************************ Loading Data **************************************
# Reading Systematic Lupus Erythematosus(SLE) data
sleData <- readRDS("240209_lupusDataSCE_raw_CM.rds")
# Assign named assay "X" to counts: For easy manupulation
counts(sleData) <- assay(sleData, "X")
dim(counts(sleData)) # 12869 307429

# Aggregation of single-cell to pseudobulk data
#************************ pseudo-bulk samples*********************************
# Creating pseudo-bulk samples: Aggregate by ind_cov_batch_cov
sle_summed <- aggregateAcrossCells(sleData, 
                                   id=colData(sleData)[,c("ind_cov_batch_cov")])
dim(counts(sle_summed)) # 12869   355
dim(colData(sle_summed)) # 355 by 15       -> 355 samples
dim(rowData(sle_summed)) # 12869  by   2   -> 12869 genes

# Exploration of the reduced Dimension slot in sle_summed:
dim(reducedDim(sle_summed, "X_pca")) # 355  50
dim(reducedDim(sle_summed, "X_umap")) # 355   2

sleData2 <- sle_summed
colData(sleData2)$mean_counts <- colMeans(assay(sleData2, "counts"))
colData(sleData2)$total_counts <- colSums(assay(sleData2, "counts"))
colData(sleData2)

#*********************** Exploration of the gene expression data*************
#  Distribution of total counts per cell is for each of our batches
ggcells(sleData2, aes(x = batch_cov, y = total_counts)) +
  geom_violin(fill = 'brown') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle("Distribution of Total Counts per Cell Across Batches")


# Distribution of total counts per cell across SLE status
ggcells(sleData2, aes(x = SLE_status, y = total_counts)) + 
  geom_violin(fill = 'orange') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle("Distribution of Total Counts per Cell  across SLE status")


# PCA plots
pcaplot1 <-plotReducedDim(sle_summed, "X_pca", colour_by = "SLE_status")
p1 <- pcaplot1 + labs(x = "PCA 1", y = "PCA 2", title = "PCA Plot by SLE Status") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

pcaplot2 <-plotReducedDim(sle_summed, "X_pca", colour_by = "Status")
p2 <- pcaplot2 + labs(x = "PCA 1", y = "PCA 2", title = "PCA Plot by Status") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

pcaplot3 <-plotReducedDim(sle_summed, "X_pca", colour_by = "Processing_Cohort")
p3 <- pcaplot3 + labs(x = "PCA 1", y = "PCA 2", title = "PCA Plot by Processing Cohort") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

pcaplot4 <-plotReducedDim(sle_summed, "X_pca", colour_by = "Sex")
p4 <- pcaplot4 + labs(x = "PCA 1", y = "PCA 2", title = "PCA Plot by Sex") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)

# UMAP Plots
umapplot1 <- plotReducedDim(sle_summed, "X_umap", colour_by = "SLE_status")
u1 <- umapplot1 + labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP Plot of SLE Status") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) 

umapplot2 <- plotReducedDim(sle_summed, "X_umap", colour_by = "Status")
u2 <- umapplot2 + labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP Plot of Status") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) 

umapplot3 <- plotReducedDim(sle_summed, "X_umap", colour_by = "Processing_Cohort")
u3<- umapplot3 + labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP Plot of Processing_Cohort") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) 

umapplot3 <- plotReducedDim(sle_summed, "X_umap", colour_by = "Sex")
u4 <- umapplot3 + labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP Plot by Sex") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) 
grid.arrange(u1, u2, u3, u4, nrow = 2, ncol = 2)

# Mean versus variance plot
dgC_Matrix <- counts(sle_summed)
mean_counts <- apply(dgC_Matrix, 1, mean)         
variance_counts <- apply(dgC_Matrix, 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  labs(title = "Mean-Variance Plot",
       x = "Mean Counts", 
       y = "Variance Counts") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


#***********************  Performing the PseudobulkDGE analysis: ***********************
#------using methods from edgeR-----------#
label <- "cM" # selecting Classical Monocyte cell type
current <- sle_summed[,label==sle_summed$cg_cov]

table(current$SLE_status) # Healthy 146; SLE 209

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(current), samples=colData(current))

# Removing  combinations that have very few or lowly-sequenced cells.
# We remove combinations containing fewer than 10 cells.
discarded <- current$ncells < 10
y <- y[,!discarded]
summary(discarded) # only one pseudobulk profile removed

# Remove genes that are lowly expressed
keep <- filterByExpr(y, group=current$SLE_status) # default min total count 15
# removes genes that are not expressed above a log-CPM threshold in a minimum number of samples 
y <- y[keep,]
summary(keep) # 8281 genes remaining; 4588 genes removed

# Correcting for composition biases by computing normalization factors.
y <- calcNormFactors(y)
y$samples

# Diagnostics for a bulk RNA-seq DE analysis

# Mean-difference (MD) plot for each normalized pseudo-bulk profile (only 6 are shown)
par(mfrow=c(2,3), mar=c(2,2,2,2))
for (i in seq_len(6)) {
  plotMD(y, column=i)
}

# multi-dimensional scaling (MDS) plot
plotMDS(cpm(y, log=TRUE), 
        col=ifelse(as.factor(y$samples$SLE_status), "blue", "red"))

# Setting up a design matrix
# Case1:
design <- model.matrix(~SLE_status, y$samples)
design

# Estimating the negative binomial (NB) dispersions
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
plotBCV(y)

# Estimating the quasi-likelihood dispersions with glmQLFit()
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)
summary(fit$df.prior)
plotQLDisp(fit)

# Testing for differences in expression due to disease effect using glmQLFTest()
# DEGs are defined as those with non-zero log-fold changes 
# at a false discovery rate of 5%.
res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res)) # Significant 5283; not significant 2998.

topTags(res)

#------------------------------------------------------------------------------
# Similar analysis as above but doing it for all cell samples;
#  Looping across labels

# Removing all pseudo-bulk samples with 'insufficient' cells.
summed.filt <- sle_summed[,sle_summed$ncells >= 10]

de.results <- pseudoBulkDGE(summed.filt, 
                            label=summed.filt$cg_cov,
                            design=~SLE_status,
                            coef="SLE_statusSLE",
                            condition=summed.filt$SLE_status 
)

cur.results <- de.results[["cM"]]
cur.results[order(cur.results$PValue),]



#  Cross-label meta-analyses
is.de <- decideTestsPerLabel(de.results, threshold=0.05)
summarizeTestsPerLabel(is.de)


# Upregulated across most cell types.
up.de <- is.de > 0 & !is.na(is.de)
head(sort(rowMeans(up.de), decreasing=TRUE), 10)

## Downregulated across cell types.
down.de <- is.de < 0 & !is.na(is.de)
head(sort(rowMeans(down.de), decreasing=TRUE), 10)

#  To identify label-specific DE; using pseudoBulkSpecific
de.specific <- pseudoBulkSpecific(summed.filt,
                                  label=summed.filt$cg_cov,
                                  design=~SLE_status,
                                  coef="SLE_statusSLE",
                                  condition=summed.filt$SLE_status
)

cur.specific <- de.specific[["cM"]]
cur.specific <- cur.specific[order(cur.specific$PValue),]
cur.specific
#-------------------------------------------------------------------------------

# Case 2: # Setting up a design matrix
design <- model.matrix(~factor(batch_cov)+ Age+Sex+pop_cov+SLE_status, y$samples) 
design

# Estimating the negative binomial (NB) dispersions
y <- estimateDisp(y, design)


# Estimating the quasi-likelihood dispersions with glmQLFit()
fit <- glmQLFit(y, design, robust=TRUE)

# Testing for differences in expression due to disease effect 
res2 <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res2)) # Significant 5283; not significant 2998.

topTags(res2)

# Plots of top 2 DE genes
sizeFactors(summed.filt) <- NULL
g1 <- plotExpression(logNormCounts(summed.filt),
                     features="CHMP5",
                     x="SLE_status", colour_by="SLE_status")

g2 <- plotExpression(logNormCounts(summed.filt),
                     features="MYL12A",
                     x="SLE_status", colour_by="SLE_status")

grid.arrange(g1, g2, nrow = 1, ncol = 2,
             top = "Expression of CHMP5 and MYL12A Across SLE Status")


#***********************  Independence assumption!!!  **************************
# NOTE: For some patients, multiple samples were taken, So no independence......
# -----------selecting one sample per patient from the aggregated data---------#
# --------- Select the samples with the highest total count --------------------# 

# Add a column 'total_counts' to the'metadata', representing total counts for each cell
sle_Data_count <- sleData
colData(sle_Data_count)$total_counts <- colSums(counts(sle_Data_count))

metadata_count <- as.data.frame(colData(sle_Data_count))

# Summarize to find the ind_cov_batch_cov with the highest total count
ind_batch_summary <- metadata_count %>%
  group_by(ind_cov, ind_cov_batch_cov) %>%
  summarise(total_counts_sum = sum(total_counts), .groups = 'drop') %>%
  arrange(ind_cov, desc(total_counts_sum)) %>%
  group_by(ind_cov) %>%
  filter(row_number() == 1) %>%
  ungroup()

dim(ind_batch_summary) # 261 by 3

# Extract the ind_cov_batch_cov identifiers to retain
ind_batch_ids_to_retain <- as.character(ind_batch_summary$ind_cov_batch_cov)

#  subset the aggregated data to contain one sample per patient
sle_summed_indep <-  sle_summed[,  sle_summed$ind_cov_batch_cov %in% ind_batch_ids_to_retain]

dim(sle_summed_indep) # 12869   261

#****Performing the PseudobulkDGE analysis: now with one sample per patient*****
#----------- using methods from edgeR --------------#
label <- "cM" # selecting Classical Monocyte cell type
sle_summed2 <- sle_summed_indep
current <- sle_summed2[,label==sle_summed2$cg_cov]

table(current$SLE_status) # Healthy 99; SLE 162

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(current), samples=colData(current))

# Removing  combinations that have very few or lowly-sequenced cells.
# We remove combinations containing fewer than 10 cells.
discarded <- current$ncells < 10
y <- y[,!discarded]
summary(discarded) # No pseudobulk profile removed => all profiles have 10 cells+

# Remove genes that are lowly expressed
keep <- filterByExpr(y, group=current$SLE_status) # default min total count 15
# removes genes that are not expressed above a log-CPM threshold in a minimum number of samples 
y <- y[keep,]
summary(keep) # 8584 genes remaining; 4285 genes removed

# Correcting for composition biases by computing normalization factors.
y <- calcNormFactors(y)
y$samples

# Diagnostics for a bulk RNA-seq DE analysis
# Mean-difference (MD) plot for each normalized pseudo-bulk profile (only 6 are shown)
par(mfrow=c(2,3), mar=c(2,2,2,2))
for (i in seq_len(6)) {
  plotMD(y, column=i)
}

# multi-dimensional scaling (MDS) plot
plotMDS(cpm(y, log=TRUE), 
        col=ifelse(as.factor(y$samples$SLE_status), "blue", "red"))

# Setting up a design matrix
# Case1:
design <- model.matrix(~SLE_status, y$samples)
design

# Estimating the negative binomial (NB) dispersions
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
plotBCV(y)

# Estimating the quasi-likelihood dispersions with glmQLFit()
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)
summary(fit$df.prior)
plotQLDisp(fit)

# Testing for differences in expression due to disease effect using glmQLFTest()
# DEGs are defined as those with non-zero log-fold changes 
# at a false discovery rate of 5%.
res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res,  adjust.method="BH", p.value=0.05, lfc=0))# Significant 4950; not significant 3634.

topTags(res)

#------------------------------------------------------------------------------
# Similar analysis as above but doing it for all cell samples;
#  Looping across labels

# Removing all pseudo-bulk samples with 'insufficient' cells.
summed.filt <- sle_summed2[,sle_summed2$ncells >= 10]

de.results <- pseudoBulkDGE(summed.filt, 
                            label=summed.filt$cg_cov,
                            design=~SLE_status,
                            coef="SLE_statusSLE",
                            condition=summed.filt$SLE_status 
)

cur.results <- de.results[["cM"]]
cur.results[order(cur.results$PValue),]



#  Cross-label meta-analyses
is.de <- decideTestsPerLabel(de.results, threshold=0.05)
summarizeTestsPerLabel(is.de)


# Upregulated across most cell types.
up.de <- is.de > 0 & !is.na(is.de)
head(sort(rowMeans(up.de), decreasing=TRUE), 10)

## Downregulated across cell types.
down.de <- is.de < 0 & !is.na(is.de)
head(sort(rowMeans(down.de), decreasing=TRUE), 10)

#  To identify label-specific DE; using pseudoBulkSpecific
de.specific <- pseudoBulkSpecific(summed.filt,
                                  label=summed.filt$cg_cov,
                                  design=~SLE_status,
                                  coef="SLE_statusSLE",
                                  condition=summed.filt$SLE_status
)

cur.specific <- de.specific[["cM"]]
cur.specific <- cur.specific[order(cur.specific$PValue),]
cur.specific
#-------------------------------------------------------------------------------

# Case 2: # Setting up a design matrix
design <- model.matrix(~factor(batch_cov)+ Age+Sex+pop_cov+SLE_status, y$samples) 
design

# Estimating the negative binomial (NB) dispersions
y <- estimateDisp(y, design)

# Estimating the quasi-likelihood dispersions with glmQLFit()
fit <- glmQLFit(y, design, robust=TRUE)

# Testing for differences in expression due to disease effect 
res2 <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res2)) # Significant 2195; not significant 6489.

topTags(res2)

# Plots of top 2 DE genes
sizeFactors(summed.filt) <- NULL
g1 <- plotExpression(logNormCounts(summed.filt),
                     features="MIR24-2",
                     x="SLE_status", colour_by="SLE_status")

g2 <- plotExpression(logNormCounts(summed.filt),
                     features="EPSTI1",
                     x="SLE_status", colour_by="SLE_status")

grid.arrange(g1, g2, nrow = 1, ncol = 2,
             top = "Expression of MIR24-2 and EPSTI1 Across SLE Status")


#******************** Fitting NB Models using glmmTMB ***********************
sce <- sle_summed_indep

# Remove genes that are lowly expressed
keep <- filterByExpr(sce, group=sce$SLE_status) # default min total count 15
# removes genes that are not expressed above a log-CPM threshold in a minimum number of samples 
sce <- sce[keep,]
summary(keep) # Remained 8584 genes

# Calculating normalised library sizes (effective library sizes)
sce_norm <- calcNormFactors(sce) #or use normLibSizes() function which is the same
NormLibSize <- sce_norm$samples$lib.size*sce_norm$samples$norm.factors

# Case1: gene expression ~ disease
# Initialize a list to store model results for all genes
models.list <- list()

# Initialize a vector to flag genes with fitting issues
problematic_genes <- character()

# Create a dataframe for glmmTMB fitting
lupus <- data.frame(
  Expression = NA,
  Disease = sce$SLE_status,
  LogLibSize = log(NormLibSize)
)

# Loop to fit models and catch any fitting problems
for (i in 1: nrow(sce)) {
  
  gene_name <- rownames(sce)[i]
  # print(paste("Processing gene:", gene_name))
  
  # Attempt to fit the model, catch warnings and errors
  fit_result <- tryCatch({
    warning_issue <- FALSE
    warning_message <- ""
    
    # Catch warnings specifically
    model <- withCallingHandlers({
      lupus$Expression <- as.numeric(counts(sce)[i, ])
      glmmTMB(Expression ~ Disease + offset(LogLibSize), data = lupus, family = nbinom2)
    }, warning = function(w) {
      warning_issue <- TRUE
      warning_message <- w$message
      invokeRestart("muffleWarning")
    })
    
    list(model = model, issue = warning_issue, message = warning_message)
  }, error = function(e) {
    # Handle errors
    list(model = NULL, issue = TRUE, message = e$message)
  })
  
  # Store the model result
  models.list[[gene_name]] <- fit_result
  
  # If there was an issue (warning or error), flag the gene
  if (fit_result$issue) {
    problematic_genes <- c(problematic_genes, gene_name)
    message_type <- ifelse(fit_result$message == "", "Error", "Warning")
    print(paste(message_type, "with gene:", gene_name, "Message:", fit_result$message))
  }
}

# Reviewing problematic_genes to see which had issues
length(problematic_genes) # How many genes had issues
problematic_genes # Which genes had issues

# Save the list of full model objects 
saveRDS(models.list, "models_list.rds")

# Load the models list from the RDS file
mlist <- readRDS("models_list.rds")


# Initialize an empty dataframe for storing Effect size, standard error and p-values
resultsNB_df <- data.frame(
  Gene = character(),
  effect_size = numeric(),
  std_error = numeric(),
  PValue = numeric() 
)

# Loop through each model in the list
for (gene in names(mlist)) {
  # Extract the model from the list
  model.fit <- mlist[[gene]]
  
  # Extract out effect size, standard error and pvalue for each gene
  effect_size <- summary(model.fit$model)$coefficients$cond[2, 1]
  std_error <- summary(model.fit$model)$coefficients$cond[2, 2]
  pvalue  <- summary(model.fit$model)$coefficients$cond[2, 4]
  
  # Append the gene name, effect size, standard error and its p-value to the dataframe
  resultsNB_df <- rbind(resultsNB_df, 
                      data.frame(Gene = gene, EffectSize = effect_size, 
                                 StdError = std_error, PValue = pvalue))
}

dim(resultsNB_df)

# View results, raw p-values
par(mar=c(4, 4, 2, 2)) 
hist(resultsNB_df$PValue, main="raw p-values", xlab="pvalue")

# Multiplicity correction using BH
resultsNB_df$p.adj<-p.adjust(resultsNB_df$PValue, method = "BH")
sum(resultsNB_df$p.adj>0.05)

# Sorting the dataframe by adjusted p-values in descending order
resultsNB_df_sorted <- resultsNB_df[order(resultsNB_df$p.adj), ]

# Viewing the top rows of the sorted dataframe
head(resultsNB_df_sorted, 10)

# Scatterplot of normalized expression of top 20 most significant genes
top20_sig_genes <- resultsNB_df %>% dplyr::arrange(p.adj) %>%
  dplyr::pull(Gene) %>% head(n=20)

# Obtaine the Log transformed the counts
sle_transform <- logNormCounts(sle_summed_indep)

# Extract the log-normalised counts for the top 20 DE genes
logcounts_data <- assay(sle_transform, "logcounts")
top20_data <- logcounts_data[rownames(logcounts_data) %in% top20_sig_genes, ]
dim(top20_data)

# Create a data frame suitable for plotting
plot_data <- as.data.frame(t(top20_data))
plot_data$Disease <- colData(sle_transform)$SLE_status
plot_data_long <- reshape2::melt(plot_data, id.vars = "Disease")


# Create the plot with ggplot2
DE20p <- ggplot(plot_data_long, aes(x = variable, y = value, color = Disease)) +
  geom_jitter(width = 0.2, size = 1) +
  labs(x = "Genes", y = "log2 Normalized Counts", 
       title = "Top 20 Significant Differentially Expressed Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
print(DE20p)

## Volcano plot

# Compute log2 fc
resultsNB_df$l2fc <- log(exp(resultsNB_df$EffectSize), base=2)

# Obtain logical vector where TRUE values denote padj values < 0.05 
# and fold change > 1.5 in either direction
resultsNB_thres <- resultsNB_df %>% 
  mutate(threshold = p.adj < 0.05 & abs(l2fc) >= 0.58)

ggplot(resultsNB_thres) +
    geom_point(aes(x = l2fc, y = -log(p.adj, base=10), colour = threshold)) +
    ggtitle("Volcano plot of SLE patients relative to Healthy controls") +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    scale_y_continuous(limits = c(0,50)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))

# Case2: gene expression ~ disease + covariates

# Initialize a list to store model results for all genes
models.list2 <- list()

# Initialize a vector to flag genes with fitting issues
problematic_genes2 <- character()


# Create a dataframe for glmmTMB fitting
lupus <- data.frame(
  Expression = NA,
  Disease = sce$SLE_status,
  batch_cov = factor(sce$batch_cov),
  Age = as.numeric(sce$Age),
  Sex = sce$Sex,
  PopCov = sce$pop_cov,
  LogLibSize = log(NormLibSize)
)


# Loop to fit models and catch any fitting problems
for (i in 1: nrow(sce)) {
  
  gene_name <- rownames(sce)[i]
  
  # Attempt to fit the model, catch warnings and errors
  fit_result <- tryCatch({
    warning_issue <- FALSE
    warning_message <- ""
    
    # Catch warnings specifically
    model <- withCallingHandlers({
      lupus$Expression <- as.numeric(counts(sce)[i, ])
      
      glmmTMB(Expression ~ Disease + batch_cov + Age + Sex +
                PopCov + offset(LogLibSize), data = lupus, family = nbinom2)
    }, warning = function(w) {
      warning_issue <- TRUE
      warning_message <- w$message
      invokeRestart("muffleWarning")
    })
    
    list(model = model, issue = warning_issue, message = warning_message)
  }, error = function(e) {
    # Handle errors
    list(model = NULL, issue = TRUE, message = e$message)
  })
  
  # Store the model result
  models.list2[[gene_name]] <- fit_result
  
  # If there was an issue (warning or error), flag the gene
  if (fit_result$issue) {
    problematic_genes <- c(problematic_genes, gene_name)
    message_type <- ifelse(fit_result$message == "", "Error", "Warning")
    print(paste(message_type, "with gene:", gene_name, "Message:", fit_result$message))
  }
}

# Error in .sql_connect_RO(.sql_dbfile(bfc)) : DB file 'local/BiocFileCache.sqlite' not found

# Reviewing problematic_genes to see which had issues
length(problematic_genes2) # How many genes had issues
problematic_genes2 # Which genes had issues

# Save the list of full model objects 
saveRDS(models.list2, "models_list2.rds")

# Load the models list from the RDS file
mlist_all <- readRDS("models_list2.rds")


# Initialize an empty dataframe for storing Effect size, standard error and p-values
resultsNB_df2 <- data.frame(
  Gene = character(),
  effect_size = numeric(),
  std_error = numeric(),
  PValue = numeric() 
)

# Loop through each model in the list
for (gene in names(mlist_all)) {
  # Extract the model from the list
  model.fit <- mlist_all[[gene]]
  
  # Extract out effect size, standard error and pvalue for each gene
  effect_size <- summary(model.fit$model)$coefficients$cond[2, 1]
  std_error <- summary(model.fit$model)$coefficients$cond[2, 2]
  pvalue  <- summary(model.fit$model)$coefficients$cond[2, 4]
  
  # Append the gene name, effect size, standard error and its p-value to the dataframe
  resultsNB_df2 <- rbind(resultsNB_df2, 
                        data.frame(Gene = gene, EffectSize = effect_size, 
                                   StdError = std_error, PValue = pvalue))
}


head(resultsNB_df2)
rows_with_na <- apply(resultsNB_df2, 1, function(x) any(is.na(x)))

# Selecting rows with NA values
df_with_na <- resultsNB_df2[!complete.cases(resultsNB_df2), ]

# Genes with NAs ??? Ignore them.
print(df_with_na)

# View results, raw p-values
par(mar=c(4, 4, 2, 2)) 
hist(resultsNB_df2$PValue, main="raw p-values", xlab="pvalue")

# Multiplicity correction using BH
resultsNB_df2$p.adj<-p.adjust(resultsNB_df2$p.adj, method = "BH")
sum(resultsNB_df2$p.adj, na.rm = TRUE)


# Sorting the dataframe by adjusted p-values in descending order
resultsNB_df2_sorted <- resultsNB_df2[order(resultsNB_df2$p.adj), ]

# Viewing the top rows of the sorted dataframe
head(resultsNB_df2_sorted, 10)


# Scatterplot of normalized expression of top 20 most significant genes
top20_sig_genes <- resultsNB_df2 %>% dplyr::arrange(p.adj) %>%
  dplyr::pull(Gene) %>% head(n=20)

# Obtaine the Log transformed the counts
sle_transform <- logNormCounts(sle_summed_indep)

# Extract the log-normalised counts for the top 20 DE genes
logcounts_data <- assay(sle_transform, "logcounts")
top20_data <- logcounts_data[rownames(logcounts_data) %in% top20_sig_genes, ]
dim(top20_data)

# Create a data frame suitable for plotting
plot_data <- as.data.frame(t(top20_data))
plot_data$Disease <- colData(sle_transform)$SLE_status
plot_data_long <- reshape2::melt(plot_data, id.vars = "Disease")


# Create the plot with ggplot2
DE20p <- ggplot(plot_data_long, aes(x = variable, y = value, color = Disease)) +
  geom_jitter(width = 0.2, size = 1) +
  labs(x = "Genes", y = "log2 Normalized Counts", 
       title = "Top 20 Significant Differentially Expressed Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
print(DE20p)

## Volcano plot

# Compute log2 fc
resultsNB_df2$l2fc <- log(exp(resultsNB_df2$EffectSize), base=2)

# Obtain logical vector where TRUE values denote padj values < 0.05 
# and fold change > 1.5 in either direction
resultsNB_thres <- resultsNB_df2 %>% 
  mutate(threshold = p.adj < 0.05 & abs(l2fc) >= 0.58)

ggplot(resultsNB_thres) +
  geom_point(aes(x = l2fc, y = -log(p.adj, base=10), colour = threshold)) +
  ggtitle("Volcano plot of SLE patients relative to Healthy controls") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

#******************** Incoperating Causal Inference ***********************
# Still reading the literature......... will start here





