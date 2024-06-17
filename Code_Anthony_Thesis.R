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
library(muscat)
library(purrr)
library(scran)
library(gridExtra)
library(glmmTMB)
library(boot)
library(ggcells)
library(VennDiagram)
library(boot.pval)
library(precmed)
library(UpSetR)
library(patchwork)

## FINAL CODE HERE ##
#************************ Loading Data **************************************
# Reading Systematic Lupus Erythematosus(SLE) data set
sleData <- readRDS("240209_lupusDataSCE_raw_CM.rds")
# Assign named assay "X" to counts: For easy manupulation
counts(sleData) <- assay(sleData, "X")

#************************ pseudo-bulk samples*********************************
# Creating pseudo-bulk samples (from cell level counts): Aggregate by ind_cov_batch_cov
sle_summed <- aggregateAcrossCells(sleData, 
                                   id=colData(sleData)[,c("ind_cov_batch_cov")])
# dim(counts(sle_summed)) # 12869 genes , 355 samples

#***********************  Independence assumption!!!  **************************
# NOTE: For some patients, multiple samples were taken, So no independence......
# -----------selecting one sample per patient from the aggregated data---------#
# --------- Select the samples with the highest total count --------------------# 

# Add a column 'total_counts' to the 'metadata', representing total counts for each cell
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

# Extract the ind_cov_batch_cov identifiers to retain
ind_batch_ids_to_retain <- as.character(ind_batch_summary$ind_cov_batch_cov)

#  Subset the aggregated data to contain one sample per patient
sle_summed_indep <-  sle_summed[, sle_summed$ind_cov_batch_cov %in% 
                                  ind_batch_ids_to_retain]

# dim(sle_summed_indep) # 12869 genes,   261 unique samples

# Delete samples from from African American and Hispanic because there are few
# data points. Only 3 from African American and 2 from Hispanic.
sle_summed_indep <- sle_summed_indep[, colData(sle_summed_indep)$pop_cov %in%
                                       c("Asian", "European")]

# Dropping unused levels in all factor columns in colData
colData(sle_summed_indep) <- DataFrame(lapply(colData(sle_summed_indep), function(x) {
  if (is.factor(x)) droplevels(x) else x
}))
# dim(sle_summed_indep) # 12869   256; NOTE: 5 samples were removed

#******************    Exploration of expression data     *************
# Mean versus variance plot
dgC_Matrix <- counts(sle_summed_indep)
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

# Un Supervised Exploratory Data Analysis
# PCA plot
pcaplot1 <-plotReducedDim(sle_summed_indep, "X_pca", colour_by = "SLE_status")
p1 <- pcaplot1 + labs(x = "PCA 1", y = "PCA 2", title = "PCA Plot by SLE Status") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

# UMAP plot
umapplot1 <- plotReducedDim(sle_summed_indep, "X_umap", colour_by = "SLE_status")
p2 <- umapplot1 + labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP Plot of SLE Status") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) 

# Combined plot
grid.arrange(p1, p2, nrow = 1, ncol = 2)

#******************** Fitting NB Models using glmmTMB ***********************
# Assign the data object 'sle_summed_indep' to the variable 'sle'to work with.
sce <- sle_summed_indep

# FILTERING
# Remove genes that are lowly expressed
keep <- filterByExpr(sce, group=sce$SLE_status) # default min total count 15
# removes genes that are not expressed above a log-CPM threshold in a minimum number of samples 
sce <- sce[keep,]; summary(keep) # Remained 8578 genes

# Calculating normalised library sizes (effective library sizes)
sce_norm <- calcNormFactors(sce, method = "TMM") # or use normLibSizes() function which is the same
NormLibSize <- sce_norm$samples$lib.size*sce_norm$samples$norm.factors

# Naive Case: gene expression ~ disease
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
      
      # Model
      glmmTMB(Expression ~ Disease + offset(LogLibSize), data = lupus, family = nbinom2(link="log"))
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


# Multiplicity correction using BH
resultsNB_df$p.adj <- p.adjust(resultsNB_df$PValue, method = "BH")
sum(resultsNB_df$p.adj > 0.05) # significant genes

# Sorting the dataframe by adjusted p-values in descending order
resultsNB_df_sorted <- resultsNB_df[order(resultsNB_df$p.adj), ]

# Viewing the top rows of the sorted dataframe
head(resultsNB_df_sorted, 10)

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

# Conditional: gene expression ~ disease + covariates
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
      
      # Model
      glmmTMB(Expression ~ Disease + batch_cov + Age + Sex +
                PopCov + offset(LogLibSize), data = lupus, family = nbinom2(link="log"))
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

# Multiplicity correction using BH
resultsNB_df2$p.adj<-p.adjust(resultsNB_df2$p.adj, method = "BH")
sum(resultsNB_df2$p.adj < 0.05, na.rm = TRUE)

# Sorting the dataframe by adjusted p-values in descending order
resultsNB_df2_sorted <- resultsNB_df2[order(resultsNB_df2$p.adj), ]

# Viewing the top rows of the sorted dataframe
head(resultsNB_df2_sorted, 10)

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

############### Venn diagram Naive case Vs Conditional case ##############
# case 1: Naive; case2 conditional

resultsNB_dfSig <- subset(resultsNB_df, p.adj < 0.05)
case1 <- resultsNB_dfSig$Gene

resultsNB_df2Sig <- subset(resultsNB_df2, p.adj < 0.05)
case2 <- resultsNB_df2Sig$Gene

# Arguments for a pairwise (two-sets) venn-diagram are sizes for set1, set2 and overlap (intersect)
venn.plot <- draw.pairwise.venn(length(case2),
                                length(case1),
                                # Calculate the intersection of the two sets
                                length(intersect(case2, case1) ),
                                category = c("Conditional Analysis", "Naive Analysis"), scaled = F,
                                fill = c("lightblue", "lightpink"), alpha = rep(0.5, 2),
                                cat.pos = c(0, 0))

# ploting 
grid.draw(venn.plot)

# -----------------  Scatter Plot comparing effect sizes ----------------------
# Filtering df1 to only include genes that are in df2
filtered_df1 <- resultsNB_df %>% filter(Gene %in% resultsNB_df2$Gene)

data_to_plot <- data.frame(Gene = filtered_df1$Gene, 
                           naive_ES = filtered_df1$EffectSize, 
                           conditional_ES = resultsNB_df2$EffectSize)

# Create the scatter plot with colored points and y = x line
ggplot(data_to_plot, aes(x = naive_ES, y = conditional_ES)) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  ggtitle("Comparison of Effect Size Estimates") +
  xlab("Naive Effect Sizes") +
  ylab("Conditional Effect Sizes") +
  theme_bw()

# -----------------  Scatter Plot comparing p-values ---------------------------
data <- data.frame(Gene = filtered_df1$Gene, 
                   Pval1 = filtered_df1$p.adj, 
                   Pval2 = resultsNB_df2$p.adj)

ggplot(data, aes(x = -log(Pval1, base = 10), y = -log(Pval2, base = 10))) +
  geom_point(alpha = 0.5, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Adding y=x line
  scale_x_log10(labels = scales::label_comma()) +  # Log scale for X-axis
  scale_y_log10(labels = scales::label_comma()) +  # Log scale for Y-axis
  ggtitle("Comparison of p-values") +
  xlab("- Log of Adj.p-values (Naive Case)") +
  ylab("- Log of Adj.p-values (Conditional Case)") +
  theme_bw()

######################## Poisson Versus Negative Binomial###############
sce <- sle_summed_indep

# Remove genes that are lowly expressed
keep <- filterByExpr(sce, group=sce$SLE_status) # default min total count 15
# removes genes that are not expressed above a log-CPM threshold in a minimum number of samples 
sce <- sce[keep,]
summary(keep) # Remained 8578 genes

# Sample 100 genes from the SingleCellExperiment object
set.seed(123)  # Set seed for reproducibility
sampled_gene_indices <- sample(nrow(sce), 100)

# Subset the SingleCellExperiment object to include only the sampled genes
sce_sampled_genes <- sce[sampled_gene_indices, ]


# Calculating normalised library sizes (effective library sizes)
sce_norm <- calcNormFactors(sce_sampled_genes, method = "TMM") #or use normLibSizes() function which is the same
NormLibSize <- sce_norm$samples$lib.size*sce_norm$samples$norm.factors

# Initialize a list to store model results for all genes

models.list2Poisson <- list()

# Initialize a vector to flag genes with fitting issues
problematic_genes2 <- character()

# Create a dataframe for glmmTMB fitting
lupus <- data.frame(
  Expression = NA,
  Disease = sce_sampled_genes$SLE_status,
  batch_cov = factor(sce_sampled_genes$batch_cov),
  Age = as.numeric(sce_sampled_genes$Age),
  Sex = sce_sampled_genes$Sex,
  PopCov = sce_sampled_genes$pop_cov,
  LogLibSize = log(NormLibSize)
)

# POISSON MODELS

# Loop to fit models and catch any fitting problems
for (i in 1: nrow(sce_sampled_genes)) {
  
  gene_name <- rownames(sce_sampled_genes)[i]
  
  # Attempt to fit the model, catch warnings and errors
  fit_result <- tryCatch({
    warning_issue <- FALSE
    warning_message <- ""
    
    # Catch warnings specifically
    model <- withCallingHandlers({
      lupus$Expression <- as.numeric(counts(sce_sampled_genes)[i, ])
      
      # Model
      glmmTMB(Expression ~ Disease + batch_cov + Age + Sex +
                PopCov + offset(LogLibSize), data = lupus,  
              family = poisson(link="log") )
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
  
  models.list2Poisson[[gene_name]] <- fit_result
  
  # If there was an issue (warning or error), flag the gene
  if (fit_result$issue) {
    problematic_genes <- c(problematic_genes, gene_name)
    message_type <- ifelse(fit_result$message == "", "Error", "Warning")
    print(paste(message_type, "with gene:", gene_name, "Message:", fit_result$message))
  }
}
# Save the list of full model objects 
saveRDS(models.list2Poisson, "models.list2Poisson.rds")

# Load the models list from the RDS file
mlist_all <- readRDS("models.list2Poisson.rds")

# Initialize an empty dataframe for storing Effect size, standard error and p-values
resultsNB_dfPoisson <- data.frame(
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
  resultsNB_dfPoisson <- rbind(resultsNB_dfPoisson, 
                               data.frame(Gene = gene, EffectSize = effect_size, 
                                          StdError = std_error, PValue = pvalue))
}

# NEGATIVE BINOMIAL MODELS
# Initialize a list to store model results for all genes

models.list2NB <- list()

# Initialize a vector to flag genes with fitting issues
problematic_genes2 <- character()

# Create a dataframe for glmmTMB fitting
lupus <- data.frame(
  Expression = NA,
  Disease = sce_sampled_genes$SLE_status,
  batch_cov = factor(sce_sampled_genes$batch_cov),
  Age = as.numeric(sce_sampled_genes$Age),
  Sex = sce_sampled_genes$Sex,
  PopCov = sce_sampled_genes$pop_cov,
  LogLibSize = log(NormLibSize)
)

# Loop to fit models and catch any fitting problems
for (i in 1: nrow(sce_sampled_genes)) {
  
  gene_name <- rownames(sce_sampled_genes)[i]
  
  # Attempt to fit the model, catch warnings and errors
  fit_result <- tryCatch({
    warning_issue <- FALSE
    warning_message <- ""
    
    # Catch warnings specifically
    model <- withCallingHandlers({
      lupus$Expression <- as.numeric(counts(sce_sampled_genes)[i, ])
      
      # Model
      glmmTMB(Expression ~ Disease + batch_cov + Age + Sex +
                PopCov + offset(LogLibSize), data = lupus,  
              family = nbinom2(link="log") )
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
  
  models.list2NB[[gene_name]] <- fit_result
  
  # If there was an issue (warning or error), flag the gene
  if (fit_result$issue) {
    problematic_genes <- c(problematic_genes, gene_name)
    message_type <- ifelse(fit_result$message == "", "Error", "Warning")
    print(paste(message_type, "with gene:", gene_name, "Message:", fit_result$message))
  }
}

# Save the list of full model objects 
saveRDS(models.list2NB, "models.list2NB.rds")

# Load the models list from the RDS file
mlist_all <- readRDS("models.list2NB.rds")

# Initialize an empty dataframe for storing Effect size, standard error and p-values
resultsNB_dfNB <- data.frame(
  Gene = character(),
  effect_size = numeric(),
  std_error = numeric(),
  PValue = numeric(),
  disp = numeric()
  )

# Loop through each model in the list
for (gene in names(mlist_all)) {
  # Extract the model from the list
  model.fit <- mlist_all[[gene]]
  
  # Extract out effect size, standard error and pvalue for each gene
  effect_size <- summary(model.fit$model)$coefficients$cond[2, 1]
  std_error <- summary(model.fit$model)$coefficients$cond[2, 2]
  pvalue  <- summary(model.fit$model)$coefficients$cond[2, 4]
  disp  <- summary(model.fit$model)$sigma
  
  # Append the gene name, effect size, standard error and its p-value to the dataframe
  resultsNB_dfNB <- rbind(resultsNB_dfNB, 
                          data.frame(Gene = gene, EffectSize = effect_size, 
                                     StdError = std_error, PValue = pvalue, disp = disp))
}

# Number of significant genes
resultsNB_dfPoisson$p.adj<-p.adjust(resultsNB_dfPoisson$PValue, method = "BH")
sum(resultsNB_dfPoisson$PValue < 0.05)# 49 significant genes

resultsNB_dfNB$p.adj<-p.adjust(resultsNB_dfNB$PValue, method = "BH")
sum(resultsNB_dfNB$PValue < 0.05) # 39 significant genes

# Set up the plotting parameters to adjust the margins
par(mar = c(2, 2, 2, 2) + 0.1)

# Combine the relevant data into a single data frame for ggplot2
data_NBvsPo <- data.frame(
  EffectSize_Poisson = resultsNB_dfPoisson$EffectSize,
  EffectSize_NB = resultsNB_dfNB$EffectSize
)

# Create the scatter plot with colored points and y = x line
ggplot(data_NBvsPo, aes(x = EffectSize_Poisson, y = EffectSize_NB)) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Effect size from Poisson Model",
    y = "Effect Size from Negative Binomial Model",
    title = "Comparison of Effect Sizes Estimates"
  ) + theme_bw()

# Bootstrap standard errors

# Create a dataframe for model fitting
lupus <- data.frame(
  Expression = NA,
  Disease = sce_sampled_genes$SLE_status,
  batch_cov = factor(sce_sampled_genes$batch_cov),
  Age = as.numeric(sce_sampled_genes$Age),
  Sex = sce_sampled_genes$Sex,
  PopCov = sce_sampled_genes$pop_cov,
  LogLibSize = log(NormLibSize)
)

# POISSON MODELS
# A function to compute the effect size;
comp_statistic <- function(data, indices) {
  subdata <- data[indices, ]
  modpoisson <-  glmmTMB(Expression ~ Disease + batch_cov + Age + Sex +
                           PopCov + offset(LogLibSize), data = subdata,  
                         family=poisson(link="log") )
  effect_size <- summary(modpoisson)$coefficients$cond[2, 1]
  
  return(effect_size)
}

poissonbootSE <- data.frame(
  Gene = character(),
  EffectSize = numeric(),
  StdError = numeric())

set.seed(123)
# loop through all genes to estimate the bootstrap standard error.
for (i in 1: nrow(sce_sampled_genes)) {
  
  gene_name <- rownames(sce_sampled_genes)[i]
  lupus$Expression <- as.numeric(counts(sce_sampled_genes)[i, ])
  boot.out <-  boot::boot(lupus, comp_statistic, 200)
  
  poissonbootSE <- rbind(poissonbootSE, 
                         data.frame(Gene = gene_name,
                                    EffectSize = boot.out$t0, 
                                    StdError = sd(boot.out$t)))
}

saveRDS(poissonbootSE, "poissonbootSE.rds")

# NB MODELS
# A function to compute the effect size;

comp_statisticNB <- function(data, indices) {
  subdata <- data[indices, ]
  modNB <-  glmmTMB(Expression ~ Disease + batch_cov + Age + Sex +
                      PopCov + offset(LogLibSize), data = subdata,  
                    family=nbinom2(link="log") )
  effect_size <- summary(modNB)$coefficients$cond[2, 1]
  
  return(effect_size)
}

NBbootSE <- data.frame(
  Gene = character(),
  EffectSize = numeric(),
  StdError = numeric())

set.seed(123)
for (i in 1: nrow(sce_sampled_genes)) {
  gene_name <- rownames(sce_sampled_genes)[i]
  lupus$Expression <- as.numeric(counts(sce_sampled_genes)[i, ])
  boot.out <-  boot::boot(lupus, comp_statisticNB, 200)
  
  NBbootSE <- rbind(NBbootSE, 
                    data.frame(Gene = gene_name,
                               EffectSize = boot.out$t0, 
                               StdError = sd(boot.out$t)))
}
saveRDS(NBbootSE, "NBbootSE.rds")

# Set up the plotting parameters to adjust the margins
par(mar = c(2, 2, 2, 2) + 0.1)

# Combine the relevant data into a single data frame for ggplot2
poissonbootSE <- readRDS("poissonbootSE.rds")
poissonbootSE$tStatPoisson <- poissonbootSE$EffectSize / poissonbootSE$StdError

NBbootSE  <- readRDS("NBbootSE.rds")
NBbootSE$tStatNB <- NBbootSE$EffectSize / NBbootSE$StdError

data_tstat <- data.frame(
  tstat_Poisson = poissonbootSE$tStatPoisson,
  tstat_NB = NBbootSE$tStatNB
)

# Create the scatter plot with colored points and y = x line
ggplot(data_tstat, aes(x = tstat_Poisson, y = tstat_NB)) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Test Statistic from Poisson Model",
    y = "Test Statistic from Negative Binomial Model",
    title = "Comparison of Test Statistics"
  ) + theme_bw()

# Looking at overdispersion parameters estimated in the NB modelling
# Creating the histogram for the distribution of the overdispersion parameters
hist(resultsNB_dfNB$disp, breaks=100, main="Distribution of the dispersion parameters",
     xlab="dispersion parameter", col="blue", border="black")
 
#******************** Incoperating Causal Inference ***********************
sle_causal <- sle_summed_indep

# Filtering
# Remove genes that are lowly expressed
keep <- filterByExpr(sle_causal, group=sle_causal$SLE_status) 
sle_causal <- sle_causal[keep,]
summary(keep) # Remained 8578 genes

# Calculating normalised library sizes (effective library sizes)
normFactors <- calcNormFactors(sle_causal)
NormLibSize <- normFactors$samples$lib.size*normFactors$samples$norm.factors

#**************************** Standardization  *****************************
# Create a dataframe
lupus <- data.frame(
  Expression = NA,
  Disease = sle_causal$SLE_status,
  batch_cov = factor(sle_causal$batch_cov),
  Age = as.numeric(sle_causal$Age),
  Sex = factor(sle_causal$Sex),
  PopCov = sle_causal$pop_cov,
  LogLibSize = log(NormLibSize)
)

# Function for prediction
predict_yHat <- function(data, mHealthy, mSLE, Status) {
  
  data$yHat <- NA 
  for (curID in 1:nrow(data)) {
    curRow <- data[curID,]
    curDisease <- Status
    
    if(curDisease == "Healthy"){
      curModel <- mHealthy
      otherModel <- mSLE
    } else {
      curModel <- mSLE
      otherModel <- mHealthy
    }
    
    beta <- summary(curModel)$coefficients[,1]
    betaOther <- summary(otherModel)$coefficients[,1]
    betaSex <- ifelse(curRow$Sex == "Male", beta["SexMale"], 0)
    betaPop <- ifelse(curRow$PopCov == "European",  beta["PopCovEuropean"], 0)
    logYhatNoBatchOrOffset <- beta["(Intercept)"] + curRow$Age * beta["Age"] + betaSex + betaPop 
    
    curBatchBetaName <- paste0("batch_cov",curRow$batch_cov)
    if(curBatchBetaName %in% names(beta)){
      betaBatch <- beta[curBatchBetaName]
    } else if (curBatchBetaName %in% names(betaOther)) {
      betaBatch <-  summary(otherModel)$coefficients[,1][curBatchBetaName]
    } else {
      betaBatch <-  0
    }
    
    logYhatNoOffset <- logYhatNoBatchOrOffset + betaBatch 
    logYhat <- logYhatNoOffset + curRow$LogLibSize
    yHat <- exp(logYhat) # Predicted outcome.
    data$yHat[curID] <- yHat
  } 
  return(data$yHat)
}

# Initialize an empty dataframe for storing Effect size estimates
resultsSTD <- data.frame(
  Gene = character(),
  effect_size = numeric()
)
# Loop to fit the gene  by gene glm models 
for (i in 1: nrow(sle_causal)) {
  
  gene_name <- rownames(sle_causal)[i] 
  lupus$Expression <- as.numeric(counts(sle_causal)[i, ])
  
  dataHealthy <- lupus[lupus$Disease == "Healthy",]
  dataDiseased <- lupus[lupus$Disease == "SLE",]
  
  mHealthy <- glm(Expression ~ batch_cov + Age + Sex + PopCov + offset(LogLibSize),
                  data = dataHealthy,
                  family = poisson(link = "log"))
  
  mSLE <-  glm(Expression ~ batch_cov + Age + Sex + PopCov + offset(LogLibSize),
               data = dataDiseased,
               family = poisson(link = "log"))
  
  Gcomp  <- cbind(YHat0 = predict_yHat(lupus, mHealthy, mSLE, "Healthy"),
                  YHat1 = predict_yHat(lupus, mHealthy, mSLE, "SLE"))
  
  Gcomp <- as.data.frame(Gcomp)
  Y.1 <- Gcomp$YHat1; Y.0 <- Gcomp$YHat0 
  logEffect <- log(mean(Y.1)/mean(Y.0)) # log ratio
  
  # Append the gene name, effect size to the dataframe
  resultsSTD <- rbind(resultsSTD, data.frame(Gene = gene_name, EffectSize = logEffect))
}
saveRDS(resultsSTD, "resultsSTD.rds")


##### Alternative  ######


#### ##########################################
# Inference using non-parametric bootstrap:

# Function to compute the causal effect
Std_effectsize <- function(data) {
  dataHealthy <- data[data$Disease == "Healthy",]
  dataDiseased <- data[data$Disease == "SLE",]
  
  mHealthy <- glm(Expression ~ batch_cov + Age + Sex + PopCov + offset(LogLibSize),
                  data = dataHealthy,
                  family = poisson(link = "log"))
  
  mSLE <-  glm(Expression ~ batch_cov + Age + Sex + PopCov + offset(LogLibSize),
               data = dataDiseased,
               family = poisson(link = "log"))
  
  Gcomp  <- cbind(YHat0 = predict_yHat(data, mHealthy, mSLE, "Healthy"),
                  YHat1 = predict_yHat(data, mHealthy, mSLE, "SLE"))
  
  Gcomp <- as.data.frame(Gcomp)
  Y.1 <- Gcomp$YHat1; Y.0 <- Gcomp$YHat0
  ATE <- log(mean(Y.1)/mean(Y.0)) # log ratio
  
  return(ATE)
}

###################### Alternative ######################
library(stdReg)
stdRegRatio <- function(data){
  data$Disease <- ifelse(data$Disease == "SLE", 1, 0)
  reg <- glm(Expression ~ Disease  + Age + PopCov   + offset(LogLibSize), 
             data = data, family = poisson(link="log"))
  reg.std <- stdGlm(fit=reg, data = data, X = "Disease", x=seq(0,1))
  out_std <- summary(reg.std, contrast = "ratio", reference=0)
  log_out_std <- log(out_std$est.table["1", "Estimate"])
  
  return(log_out_std)
  
}

###########################################################

# Bootstraping for Inference
# Create a function for bootstrapping
std_boot <- function(data, R) {
  
  # Estimate treatment effect for original data
  original_estimate <- Std_effectsize(data)

  # observed test statistic
  t.observed <- Std_effectsize(data)/sd(data$Expression)
  
  # Initialize vector to store bootstrap estimates
  boot_estimates <- vector(length = R)
  t.boot_null <- vector(length = R)
  
  for (i in 1:R) {
    
    # Resample data with replacement
    indices <- sample(1:nrow(data), nrow(data), replace = TRUE)
    resampled_data <- data[indices, ]
    resampled_data_null <- resampled_data
    nSLE <- sum(resampled_data$Disease == "SLE")
    
    # Resampling under the null-- no disease effect
    resampled_data_null$Disease <- sample(rep(c("SLE", "Healthy"),
                                              times = c(nSLE, nrow(data)- nSLE)))
    
    
    # Estimate treatment effects for resampled data
    boot_estimates[i] <- Std_effectsize(resampled_data)
    t.boot_null[i] <- Std_effectsize(resampled_data_null)/sd(resampled_data_null$Expression)
  }
  
  # CI
  ub <- quantile(boot_estimates, prob = 0.975, na.rm = TRUE)
  lb <- quantile(boot_estimates, prob = 0.025, na.rm = TRUE)
  StdError <- sd(boot_estimates, na.rm = TRUE)
  # Calculate p-value
  p_value <- (1+sum(abs(t.boot_null) >= abs(t.observed))) / (1 + length(t.boot_null))
  
  return(list(EffectSize = original_estimate,
              p_value = p_value, 
              lb = lb, ub = ub, StdError = StdError))
}

# Initialize the list to store results for each gene
resultsInferenceSTD <- list()

# Set seed for reproducibility
set.seed(12345)

# Loop over genes
for (i in 1:nrow(sle_causal)) {
  # Select gene name
  gene_name <- rownames(sle_causal)[i]
  lupus$Expression <- as.numeric(counts(sle_causal)[i, ])
  
  # Perform the bootstrap for the gene: Use 200 bootsrap replicates.
  bootResults <- std_boot(lupus, R=200)
  
  # Store results in a list for this gene
  geneResults <- list(
    EffectSize = bootResults$EffectSize,
    StdError = bootResults$StdError,
    lb = bootResults$lb,
    ub = bootResults$ub,
    PValue = bootResults$p_value
  )
  
  resultsInferenceSTD[[gene_name]] <- geneResults
}

# Significant genes
# Multiplicity correction using BH
resultsInferenceSTD$p.adj <- p.adjust(resultsInferenceSTD$PValue, method = "BH")

# Significant genes after multiplicity correction
sum(resultsInferenceSTD$p.adj < 0.05, na.rm = TRUE)# 3701

#****************** Inverse Probability Weighting *******************************

# Fit a logistic regression model of exposure since we have a binary treatment 0;1
lupusIPW <- data.frame(
  Expression = NA,
  Disease = sle_causal$SLE_status,
  batch_cov = factor(sle_causal$batch_cov),
  Age = as.numeric(sle_causal$Age),
  Sex = sle_causal$Sex,
  PopCov = sle_causal$pop_cov,
  LogLibSize = log(NormLibSize)
)

# Compute weights and add them (as -- sw --)  to the data frame
lupusIPW$Disease = ifelse(lupusIPW$Disease == "SLE", 1, 0)
denom.fit <- glm(Disease ~ batch_cov + Age + Sex + PopCov, 
                 family = binomial(link="logit"), data = lupusIPW)
denom.ps <- predict(denom.fit, type = "response") # denominator propensity score
numer.fit <- glm(Disease ~ 1, family = binomial(link="logit"), data = lupusIPW)
numer.ps <- predict(numer.fit, type = "response") # numerator score

# Add IP (standardized) weights to the dataset
lupusIPW$sw <- ifelse(lupusIPW$Disease == 0, 
                      ((1-numer.ps)/(1-denom.ps)), (numer.ps/denom.ps))

# Fit a glm NB model for each gene with the weights
models.listIPW <- list()
# Loop to fit the gene  by gene glm models 
ngenes <- nrow(sle_causal)
for (i in 1: ngenes) {
  
  gene_name <- rownames(sle_causal)[i]
  
  lupusIPW$Expression <- as.numeric(counts(sle_causal)[i, ])
  
  # Fit a weighted regression model
  modIPW <- glm(Expression ~ Disease + offset(LogLibSize), data = lupusIPW,
                    weights = sw, poisson(link = "log"))
  
  models.listIPW[[gene_name]] <-  modIPW
}

# Initialize an empty dataframe for storing Effect size, standard error and p-values
resultsIPW <- data.frame(
  Gene = character(),
  ATE = numeric(),
  std_error = numeric(),
  PValue = numeric() 
)

# Loop through each model in the list
for (gene in names(models.listIPW)) {
  # Extract the model from the list
  model.fit <- models.listIPW[[gene]]
  
  # Extract out effect size, standard error and pvalue for each gene
  ATE <- summary(model.fit)$coefficients[2, 1] # Effect size on a log scale.
  std_error <- summary(model.fit)$coefficients[2, 2]
  pvalue  <- summary(model.fit)$coefficients[2, 4]
  
  # Append the gene name, effect size, standard error and its p-value to the dataframe
  resultsIPW <- rbind(resultsNB_IPW, 
                         data.frame(Gene = gene, ATE = ATE, 
                                    StdError = std_error, PValue = pvalue))
}

saveRDS(resultsIPW, "resultsIPW.rds")

#  Bootstrapping for Inference.

# Create function IPW.Statistic to compute the statistic Ratio_Effect
IPW.Effectsize <- function(data) {
  # Fit a weighted regression model
  modIPW <- glm(Expression ~ Disease + offset(LogLibSize), data = data,
                weights = sw, family = poisson(link = "log"))
  EffectRatio <- summary(modIPW)$coefficients[2, 1] # Effect ratio on a log scale
  return(EffectRatio)
}

# Create a function for bootstrapping
IPW_boot <- function(data, R) {
  # Estimate treatment effect for original data
  original_estimate <- IPW.Effectsize(data)
  t.observed <- IPW.Effectsize(data)/sd(data$Expression)
  
  # Initialize vector to store bootstrap estimates
  boot_estimates <- vector(length = R)
  t.boot_null <- vector(length = R)
  
  n = nrow(data)
  for (i in 1:R) {
    # Resample data with replacement
    indices <- sample(1:n, n, replace = TRUE)
    resampled_data <- data[indices, ]
    resampled_data_null <- resampled_data
    nSLE <- sum(resampled_data$Disease == 1)
    
    # Assign disease status randomly to resample under the null (for Pvalue)
    resampled_data_null$Disease <- sample(rep(c(1, 0), 
                                              times = c(nSLE, nrow(data)- nSLE)))
    
    boot_estimates[i] <- IPW.Effectsize(resampled_data)
    t.boot_null[i] <- IPW.Effectsize(resampled_data_null)/sd(resampled_data_null$Expression)
    
  }
  
  # calculate std error
  StdError <- sd(boot_estimates)
  
  # calculate CI
  ub <- quantile(boot_estimates, prob = 0.975)
  lb <- quantile(boot_estimates, prob = 0.025)
  
  # Calculate p-value
  p_value <- (1+sum(abs(t.boot_null) >= abs(t.observed))) /(1 + length(t.boot_null))
  
  return(list(EffectSize = original_estimate,
              p_value = p_value, 
              StdError = StdError,
              ub = ub,
              lb = lb))
}

# Initialize the list to store results for each gene
resultsInferenceIPW <- list()

# Set seed for reproducibility
set.seed(12345)

# Loop over genes
for (i in 1: nrow(sle_causal)) {
  # Select gene name
  gene_name <- rownames(sle_causal)[i]
  # lupusIPW$Expression  <- as.numeric(counts(sle_causal[rownames(sle_causal) == "CXCL10", ]))
  
  lupusIPW$Expression <- as.numeric(counts(sle_causal)[i, ])
  
  # Perform the bootstrap for the gene
  bootResults <- IPW_boot(lupusIPW, R=200)
  
  # Store results in a list for this gene
  geneResults <- list(
    EffectSize = bootResults$EffectSize,
    StdError = bootResults$StdError,
    lb = bootResults$lb,
    ub = bootResults$ub,
    PValue = bootResults$p_value)
  
  resultsInferenceIPW[[gene_name]] <- geneResults
}
saveRDS(resultsInferenceIPW, "IPW_list.rds")

# Load the models list from the RDS file
result_list <- readRDS("IPW_list.rds")

# Initialize an empty dataframe for storing Effect size, standard error and p-values
resultsIPWboot <- data.frame(
  Gene = character(),
  effect_size = numeric(),
  std_error = numeric(),
  lb = numeric(),
  up = numeric(),
  PValue = numeric() 
)

# Loop through each model in the list
for (gene in names(result_list)) {
  # Extract the model from the list
  res <- result_list[[gene]]
  
  # Extract out effect size, standard error and pvalue for each gene
  effect_size <- res$EffectSize
  std_error <- res$StdError
  lb <- res$lb
  ub <- res$ub
  pvalue  <- res$PValue
  
  # Append the gene name, effect size, standard error and its p-value to the dataframe
  resultsIPWboot <- rbind(resultsIPWboot, 
                      data.frame(Gene = gene, EffectSize = effect_size, 
                                 StdError = std_error, lb = lb, ub = ub,
                                 PValue = pvalue))
}

# Significant genes
# Multiplicity correction using BH
resultsIPWboot$p.adj <- p.adjust(resultsIPWboot$PValue, method = "BH")

# Significant genes after multiplicity correction
sum(resultsIPWboot$p.adj < 0.05, na.rm = TRUE)# 3641

# Comparison between IPW Vs Standardization

## Number of significant genes ---- All significant genes
IPWSig <- subset(resultsIPWboot, p.adj < 0.05)
STWSig <- subset(resultsInferenceSTD, p.adj < 0.05)

stdgenes <- STWSig$Gene
ipwgenes <- IPWSig$Gene

venn.plot <- draw.pairwise.venn(length(ipwgenes),
                                length(stdgenes),
                                # Calculate the intersection of the two sets
                                length( intersect(ipwgenes, stdgenes) ),
                                category = c("IPW", "Standardization"), scaled = FALSE,
                                fill = c("lightblue", "lightpink"), alpha = rep(0.5, 2),
                                cat.pos = c(0, 0))
grid.draw(venn.plot)

## Number of significant genes ---- Only top 100 significant genes
IPWSigsorted <- IPWSig[order(IPWSig$p.adj, decreasing = FALSE), ]
STWSigsorted <- STWSig[order(STWSig$p.adj, decreasing = FALSE), ]

ipw2 <- head(IPWSigsorted, 100)
std2 <- head(STWSigsorted, 100)

stdgenes100 <- std2$Gene
ipwgenes100 <- ipw2$Gene

venn.plot <- draw.pairwise.venn(length(ipwgenes100),
                                length(stdgenes100),
                                # Calculate the intersection of the two sets
                                length( intersect(ipwgenes100, stdgenes100) ),
                                category = c("IPW", "Standardization"), scaled = FALSE,
                                fill = c("lightblue", "lightpink"), alpha = rep(0.5, 2),
                                cat.pos = c(0, 0))
grid.draw(venn.plot)

## Effect size estimates
load("resultsSTD.RData")
load("resultsIPW.RData")

# Combine the relevant data into a single data frame for ggplot2
IPW_Std <- data.frame(
  EffectSize_IPW = resultsIPW$ATE,
  EffectSize_STD = resultsSTD$EffectSize
)

# Create the scatter plot with colored points and y = x line
ggplot(IPW_Std, aes(x = EffectSize_IPW, y = EffectSize_STD)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Effect sizes IPW",
       y = "Effect Sizes Standardization",
       title = "Comparison of Effect Sizes Estimates") + theme_bw()

#********** Doubly Robust Approach: atefitcount using precmed package ***********
dataLupus <- data.frame(
  Expression = NA,
  Disease = sle_causal$SLE_status,
  batch_cov = factor(sle_causal$batch_cov),
  Age = as.numeric(sle_causal$Age),
  Sex = sle_causal$Sex,
  PopCov = sle_causal$pop_cov,
  LogLibSize = log(NormLibSize)
)

dataLupus$Disease <- ifelse(dataLupus$Disease == "SLE", 1, 0)

# Function to compute the effect size estimate
ateLupus <- function(data) {
  # Use tryCatch to handle errors 
  tryCatch({
    # Fit the model using the catefitcount function from the precmed package
    fit <- precmed::atefitcount(data = data,
                                cate.model = Expression ~ batch_cov + Age + Sex + PopCov + offset(LogLibSize),
                                ps.model = Disease ~ batch_cov + Age + Sex + PopCov,
                                seed = 999, verbose = 0, n.boot = 200)
    
    # Calculate the score from the fitted model
    score <- fit$log.rate.ratio
    
    # Return the computed score
    return(score)
  }, error = function(e) {
    # Return NA if any
    return(NA)
  })
}
resultsateRatio <- data.frame(
  Gene = character(),
  effect_size = numeric(),
  se = numeric(),
  lb = numeric(),
  ub = numeric(),
  pvalue = numeric())

for (i in 1: nrow(sle_causal)) {
  
  gene_name <- rownames(sle_causal)[i] 
  
  dataLupus$Expression <- as.numeric(counts(sle_causal)[i, ])
  result <- ateLupus(dataLupus)
  
  effect_size <- result[,1]
  se <- result[, 2]
  lb <- result[, 3]
  ub <- result[, 4]
  pvalue <- result[, 5]
  resultsateRatio <- rbind(resultsateRatio, 
                           data.frame(Gene = gene_name, effect_size = effect_size,
                                      se = se, lb = lb, ub = ub, pvalue = pvalue))
}

saveRDS(resultsateRatio, "resultsateRatio.rds")
resultsateRatio <- readRDS("resultsateRatio.rds")

# Significant genes
sum(resultsateRatio$pvalue < 0.05, na.rm = TRUE) # 3379

# Multiplicity correction using BH
resultsateRatio$p.adj<-p.adjust(resultsateRatio$pvalue, method = "BH")

# Significant genes after multiplicity correction
sum(resultsateRatio$p.adj < 0.05, na.rm = TRUE)# 2609


# Sorting the dataframe by adjusted p-values in descending order
resultsateRatio_sorted <- resultsateRatio[order(resultsNB_df2$p.adj), ]

# Viewing the top rows of the sorted dataframe
head(resultsateRatio_sorted, 10)

# Comparing STD vs Double robust|||| IPW vs Double Robust

## Effect size estimates
# Combine the relevant data into a single data frame for ggplot2
IPW_Std_DR <- data.frame(
  EffectSize_IPW = resultsIPW$ATE,
  EffectSize_STD = resultsSTD$EffectSize,
  EffectSize_DR = df_combined$effect_size
)

# Create the scatter plot; IPW Vs DR
p1 <- ggplot(IPW_Std_DR, aes(x = EffectSize_IPW, y = EffectSize_DR)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "IPW Effect Sizes",
       y = "Doubly Robust Effect Sizes",
       title = "Doubly Robust Vs IPW") + theme_bw()

# Create the scatter plot; STD Vs DR
p2 <- ggplot(IPW_Std_DR, aes(x = EffectSize_STD, y = EffectSize_DR)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Standardization Effect Sizes",
       y = "Doubly Robust Effect Sizes",
       title = "Doubly Robust Vs Standardization") + theme_bw()

grid.arrange(p1, p2, nrow = 1, ncol = 2)

## Number of significant genes --- All the the genes

IPWSig <- subset(resultsIPWboot, p.adj < 0.05)
STWSig <- subset(resultsInferenceSTD, p.adj < 0.05)
DRSig <- subset(resultsateRatio, p.adj < 0.05)


# gene vectors
stdgenes <- STWSig$Gene
ipwgenes <- IPWSig$Gene
DRgenes <- DRsig$Gene

# Combine all gene vectors and create a unique list of genes
allGenes <- unique(c(stdgenes, ipwgenes, DRgenes))

# Create a binary matrix indicating presence of each gene in each method
genePresence <- data.frame(
  STD = as.integer(allGenes %in% stdgenes),
  IPW = as.integer(allGenes %in% ipwgenes),
  DR = as.integer(allGenes %in% DRgenes)
)

# Set row names to gene names for clarity
rownames(genePresence) <- allGenes

# Create the UpSet plot
upset(genePresence, sets = c("STD", "IPW", "DR"))

upset(
  genePresence,
  sets = c("STD", "IPW", "DR"),
  #main.bar.color = "#56B4E9",
  matrix.color = "#D55E00",
  sets.bar.color = "#009E73",
  number.angles = 30, # Adjusts the angle of numbers
  point.size = 3,
  line.size = 2,
  keep.order = TRUE,
  # set_size.show = TRUE
)

## Number of significant genes --- Only top 100 significant genes

# Comparison between IPW, STD and Doubly Robust

DRSigsorted <- DRSig[order(DRSig$p.adj, decreasing = FALSE), ]
DR2 <- head(DRSigsorted, 100)

stdgenes100 <- std2$Gene
ipwgenes100 <- ipw2$Gene
DRgenes100  <- DR2$Gene

# Combine all gene vectors and create a unique list of genes
allGenes <- unique(c(stdgenes100, ipwgenes100, DRgenes100))

# Create a binary matrix indicating presence of each gene in each method
genePresence <- data.frame(
  STD = as.integer(allGenes %in% stdgenes100),
  IPW = as.integer(allGenes %in% ipwgenes100),
  DR = as.integer(allGenes %in% DRgenes100)
)

# Set row names to gene names for clarity
rownames(genePresence) <- allGenes

# Create the UpSet plot
upset(genePresence, sets = c("STD", "IPW", "DR"))

upset(
  genePresence,
  sets = c("STD", "IPW", "DR"),
  #main.bar.color = "#56B4E9",
  matrix.color = "#D55E00",
  sets.bar.color = "#009E73",
  number.angles = 30, # Adjusts the angle of numbers
  point.size = 3,
  line.size = 2,
  keep.order = TRUE,
  # set_size.show = TRUE
)
#********************          Yadlowsky Methods            **************

# ------------- contrastReg ------------ #

# Select 6 Genes 
genes <- c("IFI27", "OAS1", "MT2A", "KIF22", "MT-ND4L", "EIF4B")

# Initialize a list to store plots
plots <- list()

# Loop through each gene and plot
for (i in seq_along(genes)) {
  gene <- genes[i]  
  
  # Subset the data for the current gene
  dataLupus$Expression <- as.numeric(counts(sle_causal[rownames(sle_causal) == gene, ]))
  
  # Fit the model
  fit <- precmed::catefitcount(data = dataLupus,
                               score.method = "contrastReg",
                               cate.model = Expression ~ batch_cov + Age + Sex + PopCov + offset(LogLibSize),
                               ps.model = Disease ~ batch_cov + Age + Sex + PopCov,
                               initial.predictor.method = "boosting",
                               higher.y = FALSE,
                               seed = 999, verbose = 0)
  
  # Collect scores for this gene
  gene_scores <- data.frame(gene = gene, score = fit$score.contrastReg, status = dataLupus$Disease)
  
  # Create the histogram for the current gene
  p <- ggplot(gene_scores, aes(x = score, fill = status)) +
    geom_histogram(position = "stack", bins = 50, alpha = 0.6) +
    ggtitle(gene) +
    labs(x = "log CATE scores", y = "Frequency") +
    scale_fill_manual(values = c("SLE" = "red", "Healthy" = "green")) +
    theme(plot.title = element_text(size = 5, face = "bold", hjust = 0.5))
  # Disable legend except for the last plot
  if (i != length(genes)) {
    p <- p + theme(legend.position = "none")
  }
  
  plots[[i]] <- p
}

# Combine all plots
combined_plot <- wrap_plots(plots, nrow = 2, ncol = 3) + 
  # plot_annotation(title = "Distribution of the estimated log CATE scores for selected set of Genes", theme = theme(plot.title = element_text(hjust = 0.5))) +
  guides(fill = guide_legend(title = "Status"))
print(combined_plot)

# ------------- twoReg ------------ #

# Select 6 Genes 
genes <- c("IFI27", "OAS1", "MT2A", "KIF22", "MT-ND4L", "EIF4B")

# Initialize a list to store plots
plots <- list()

# Loop through each gene and plot
for (i in seq_along(genes)) {
  gene <- genes[i]  
  
  # Subset the data for the current gene
  dataLupus$Expression <- as.numeric(counts(sle_causal[rownames(sle_causal) == gene, ]))
  
  # Fit the model
  fit <- precmed::catefitcount(data = dataLupus,
                               score.method = "twoReg",
                               cate.model = Expression ~ batch_cov + Age + Sex + PopCov + offset(LogLibSize),
                               ps.model = Disease ~ batch_cov + Age + Sex + PopCov,
                               initial.predictor.method = "boosting",
                               higher.y = FALSE,
                               seed = 999, verbose = 0)
  
  # Collect scores for this gene
  gene_scores <- data.frame(gene = gene, score = fit$score.twoReg, status = dataLupus$Disease)
  
  # Create the histogram for the current gene
  p <- ggplot(gene_scores, aes(x = score, fill = status)) +
    geom_histogram(position = "stack", bins = 50, alpha = 0.6) +
    ggtitle(gene) +
    labs(x = "log CATE scores", y = "Frequency") +
    scale_fill_manual(values = c("SLE" = "red", "Healthy" = "green")) +
    theme(plot.title = element_text(size = 5, face = "bold", hjust = 0.5))
  # Disable legend except for the last plot
  if (i != length(genes)) {
    p <- p + theme(legend.position = "none")
  }
  
  plots[[i]] <- p
}

# Combine all plots
combined_plot <- wrap_plots(plots, nrow = 2, ncol = 3) + 
  # plot_annotation(title = "Distribution of the estimated log CATE scores for selected set of Genes", theme = theme(plot.title = element_text(hjust = 0.5))) +
  guides(fill = guide_legend(title = "Status"))
print(combined_plot)

# Comparison of the Effect Estimates using density plots:

# Define the  Genes 
genes <- c("IFI27", "OAS1", "MT2A", "KIF22", "MT-ND4L", "EIF4B")

# Create a new plot window with the desired layout
par(mfrow = c(2, 3))

# Loop through each gene and plot
for (i in seq_along(genes)) {
  gene <- genes[i]
  
  # Subset the data for the current gene
  dataLupus$Expression <- as.numeric(counts(sle_causal[rownames(sle_causal) == gene, ]))
  
  # Fit the model
  output_catefit <- precmed::catefitcount(data = dataLupus,
                                          score.method = c("contrastReg", "twoReg"),
                                          cate.model = Expression ~ batch_cov + Age + Sex + PopCov + offset(LogLibSize),
                                          ps.model = Disease ~ batch_cov + Age + Sex + PopCov,
                                          initial.predictor.method = "boosting",
                                          higher.y = FALSE,
                                          seed = 999, verbose = 0)
  
  dataplot <- data.frame(score = factor(rep(c("Contrast regression", "Two regressions"), 
                                            each = length(output_catefit$score.twoReg))), 
                         value = c(output_catefit$score.contrastReg, output_catefit$score.twoReg))
  
  par(mar = c(3, 3, 2, 1)) # Adjust the margin
  
  # Plot histogram for the current gene
  assign(paste0("p", i), 
         dataplot %>% 
           ggplot(aes(x = value, fill = score)) + 
           geom_density(alpha = 0.5) +
           theme_classic() + 
           labs(x = "Estimated CATE score", y = "Density", fill = "Method") +
           ggtitle(gene))
}

# Print out the stored plots
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2, 
             top = "Contrast Regression Versus Two Regression log CATE scores for selected set of Genes")

#####################################  END #################################
