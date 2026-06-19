#!/usr/bin/env Rscript

# Usage: Rscript eval_time_pred.R --train_species {STR} --target_species {STR} [OPTIONAL: --cell_type {STR}]

# ----------------- load required packages -----------------
library(data.table)
library(ggplot2)
library(stringr)
options(bitmapType='cairo') #to address the X11 issue
dodge <- position_dodge(width = 0.9)
library(VGAM) # for kendall tau calculation
library(Rmisc)
library(dplyr)
library(optparse)

# ----------------- set variables with arguments -----------------
# Define the allowed options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "input .pth file", metavar = "character"),
  make_option(c("-p", "--params"), type = "character", default = NULL, 
              help = "string with best params lr_nlayer_ndim", metavar = "character"),
  make_option(c("-x", "--train_species"), type = "character", default = NULL, 
              help = "species used for model training", metavar = "character"),
  make_option(c("-y", "--target_species"), type = "character", default = NULL, 
              help = "species used for target pseudotime prediction", metavar = "character"),
  make_option(c("-c", "--cell_type"), type = "character", default = "", 
              help = "column name used for cell type", metavar = "character"),
  make_option(c("-a", "--age"), type = "character", default = "", 
              help = "column name used for age/time of target species", metavar = "character")
)

# Parse the incoming command line parameters
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Enforce mandatory fields
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Missing mandatory argument: --input", call. = FALSE)
}
if (is.null(opt$params)){
  print_help(opt_parser)
  stop("Missing mandatory argument: --params", call. = FALSE)
}
if (is.null(opt$train_species)){
  print_help(opt_parser)
  stop("Missing mandatory argument: --train_species", call. = FALSE)
}
if (is.null(opt$target_species)){
  print_help(opt_parser)
  stop("Missing mandatory argument: --target_species", call. = FALSE)
}

input           <- opt$input
params          <- opt$params
train_species   <- opt$train_species
target_species  <- opt$target_species
cell_type       <- opt$cell_type
age             <- opt$age


# ----------------- Functions -----------------
## load the prediction file
load_existing_prediction <- function(prefix, suffix){
  time_pred_cavebear_best <- c()
  for (learning_rate in c('0.1','0.01', '0.001', '0.0001')){
    for (embed_dim in c(50, 100, 200, 400, 800)){
      for (nlayer in c(2,3,4)){
        bear_filename <- paste0(prefix, learning_rate, '_', nlayer, '_', embed_dim, suffix)
        if (file.exists(bear_filename)){
          print(bear_filename)
          time_pred_cavebear_best <- fread(bear_filename)
        }
      }
    }
  }
  return(time_pred_cavebear_best)
}

## calculate AUROC between any pair of time points, summarized by time gap
calc_pairwise_auroc <- function(input_rank, input_pred){
  # input_rank: rank of time labels (held out in target species)
  # input_pred: time prediction
  output_mat <- c()
  ntime = max(input_rank)
  for (gap in 1:(ntime-1)){
    for (i in 1:(ntime-gap)){
      j = i+gap
      if (j <= ntime){
        # make correct pred as 1 and incorrect as 0, calculate fraction of times the prediction is correct using mean
        mat <- sign(outer(input_pred[input_rank==j], input_pred[input_rank==i], FUN = "-"))
        mat[mat<0] = 0
        result <- mean(mat)
        output_mat <- rbind(output_mat, c(gap, i, result))
      }
    }
  }
  output_mat <- as.data.table(output_mat)
  return(output_mat)
}


## compare auroc pairwise by cell_type
compare_pairwise_auroc_celltype <- function(celltype_large, time_pred_cavebear, plot_dir){
  print(plot_dir)
  auc_mat <- c()
  for (celltype_i in celltype_large){
    time_pred_celltypei <- time_pred_cavebear[cell_type==celltype_i]
    time_pred_celltypei <- time_pred_celltypei[!is.na(age)]
    time_pred_celltypei$timerank <- dense_rank(time_pred_celltypei$age)
    
    ## calculate pair-wise AUROC on each cell type. try different cell type level to test the robustness of Cavebear
    for (celltype_sub_cate in c(cell_type)){
      time_pred_celltypei$celltype_sub_cate <- time_pred_celltypei[[celltype_sub_cate]]
      for (celltype_sub_i in unique(time_pred_celltypei$celltype_sub_cate)){
        time_pred_celltypei_sub <- time_pred_celltypei[celltype_sub_cate==celltype_sub_i]
        Cavebear_mat <- calc_pairwise_auroc(time_pred_celltypei_sub$timerank, time_pred_celltypei_sub$Cavebear)
        names(Cavebear_mat) <- c('gap', 'timerank', 'Cavebear')
        Cavebear_mat$cell_type <- celltype_sub_i
        Cavebear_mat$celltype_sub_cate <- celltype_sub_cate
        auc_mat <- rbind(auc_mat, Cavebear_mat)
      }
    }
  }

  auc_mat_melt <- melt(auc_mat, id.vars = c('gap', 'cell_type', 'celltype_sub_cate'), measure.vars = c('Cavebear'), variable.name = c('method'), value.name = "precision")
  auc_mat_melt$gap <- factor(auc_mat_melt$gap, levels=unique(auc_mat_melt$gap)[order(unique(auc_mat_melt$gap))])
  auc_mat_melt_stat <- summarySE(auc_mat_melt, measurevar="precision", groupvars=c("gap", "method", "cell_type"))
  auc_mat_melt_stat$method <- factor(auc_mat_melt_stat$method, levels=c('Cavebear'))
  
  p_compare_meanline_bearonly  <- ggplot(auc_mat_melt_stat[auc_mat_melt_stat$method=='Cavebear',], aes(x=gap, y=precision, colour=method)) + facet_wrap(~cell_type) +
    geom_errorbar(aes(ymin=precision-se, ymax=precision+se), width=.1) +
    geom_line() +
    geom_point()+
    xlab('Time gap') +
    ylab('Accuracy') +
    geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
    #theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
    theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
    theme(axis.text.x = element_text(angle = 75, vjust = 0.5))
  
  ggsave(p_compare_meanline_bearonly, file=paste0(plot_dir, 'compare_pseudotime_alltime_bear_byct_barplot.png'), width=12, height=6)

}

calc_pairwise_auroc_species <- function(time_pred_cavebear, plot_dir){
  print(plot_dir)
  auc_mat <- c()

  time_pred_cavebear$timerank <- dense_rank(time_pred_cavebear$age)
  Cavebear_mat <- calc_pairwise_auroc(time_pred_cavebear$timerank, time_pred_cavebear$Cavebear)
  names(Cavebear_mat) <- c('gap', 'timerank', 'Cavebear')
  auc_mat <- rbind(auc_mat, Cavebear_mat)

  auc_mat_melt <- melt(auc_mat, id.vars = c('gap'), measure.vars = c('Cavebear'), variable.name = c('method'), value.name = "precision")
  auc_mat_melt$gap <- factor(auc_mat_melt$gap, levels=unique(auc_mat_melt$gap)[order(unique(auc_mat_melt$gap))])
  auc_mat_melt_stat <- summarySE(auc_mat_melt, measurevar="precision", groupvars=c("gap", "method"))
  auc_mat_melt_stat$method <- factor(auc_mat_melt_stat$method, levels=c('Cavebear'))
  
  p_compare_meanline_bearonly  <- ggplot(auc_mat_melt_stat[auc_mat_melt_stat$method=='Cavebear',], aes(x=gap, y=precision, colour=method)) +
    geom_errorbar(aes(ymin=precision-se, ymax=precision+se), width=.1) +
    geom_line() +
    geom_point()+
    xlab('Time gap') +
    ylab('Accuracy') +
    geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
    #theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
    theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
    theme(axis.text.x = element_text(angle = 75, vjust = 0.5))
  
  ggsave(p_compare_meanline_bearonly, file=paste0(plot_dir, 'compare_pseudotime_alltime_bear_barplot.png'), width=12, height=6)

}


## ==============================================================
## compare pairwise AUROC, trained and validated across all time points in target species
## ==============================================================
input_dir <- dirname(input)
name <- sub(pattern = "\\.h5ad$", replacement = "", x = basename(opt$input))
plot_dir <- paste0(dirname(input),'/eval/')
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

## load Cavebear prediction
prefix <- paste0(input_dir, '/cvae_pytorch_disc_best_model_', name,'_', params, '_')
suffix <- paste0('_',train_species,'_pred_time',target_species,'.txt')
time_pred_cavebear <- load_existing_prediction(prefix, suffix)
names(time_pred_cavebear)[ncol(time_pred_cavebear)] <- 'Cavebear'
if (cell_type != "") {
  time_pred_cavebear <- time_pred_cavebear %>% rename("cell_type" := !!sym(cell_type))
  cell_type <- "cell_type"  # ← add this
}
if (age != "") {
  time_pred_cavebear <- time_pred_cavebear %>% rename("age" := !!sym(age))
  age <- "age"  # ← add this
}
setDT(time_pred_cavebear)

if (cell_type != "") {
  celltype_large <- names(table(time_pred_cavebear$cell_type))[table(time_pred_cavebear$cell_type)>10000] 
  compare_pairwise_auroc_celltype(celltype_large, time_pred_cavebear, plot_dir)
} else {
  compare_pairwise_auroc_species(time_pred_cavebear, plot_dir)
}
