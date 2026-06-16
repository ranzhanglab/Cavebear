#!/usr/bin/env Rscript

# Usage: Rscript eval_time_pred.R {species} {input_}

library(data.table)
library(ggplot2)
library(stringr)
options(bitmapType='cairo') #to address the X11 issue
dodge <- position_dodge(width = 0.9)
library(VGAM) # for kendall tau calculation
library(Rmisc)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

target_species = args[1] # options: zebrafish

## calculate AUROC between any pair of time points, summarized by time gap
calc_pairwise_auroc <- function(input_rank, input_pred){
  # input_rank: rank of time labels (held out in target species)
  # input_pred: time prediction
  output_mat <- c()
  ntime = max(input_rank)
  for (gap in 1:(ntime-1)){
    boolean_vec <- c()
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


compare_pairwise_auroc <- function(celltype_large, time_pred_cavebear, monocle_input_dir, monocle_train_ver, cavebear_input_dir, plot_dir){
  print(plot_dir)
  auc_mat <- c()
  for (celltype_i in celltype_large){
    celltype_i_name <- gsub('/', '_', celltype_i)
    filename <- paste0(monocle_input_dir, 'monocle2/', 'pseudotime_', monocle_train_ver, celltype_i_name, '.txt')
    if (file.exists(filename)){
      ## load prediction, find overlap between Cavebear and Monocle 2, and remove cells with na or infinite monocle predictions
      time_pred_monocle <- fread(filename)
      time_pred_monocle_combined <- merge(time_pred_cavebear, time_pred_monocle, by.x='cell_unique', by.y = 'cell')
      time_pred_combined_celltypei <- time_pred_monocle_combined[major_trajectory==celltype_i]
      time_pred_combined_celltypei <- time_pred_combined_celltypei[!is.infinite(monocle)]
      time_pred_combined_celltypei <- time_pred_combined_celltypei[!is.na(timepoint)]
      time_pred_combined_celltypei$timerank <- dense_rank(time_pred_combined_celltypei$age)
      
      ## calculate pair-wise AUROC on each cell type. try different cell type level to test the robustness of Cavebear
      for (celltype_sub_cate in c('major_trajectory', 'predicted_cell_type', 'cell_type_broad', 'cell_type_sub')){
        time_pred_combined_celltypei$celltype_sub_cate <- time_pred_combined_celltypei[[celltype_sub_cate]]
        for (celltype_sub_i in unique(time_pred_combined_celltypei$celltype_sub_cate)){
          time_pred_combined_celltypei_sub <- time_pred_combined_celltypei[celltype_sub_cate==celltype_sub_i]
          Cavebear_mat <- calc_pairwise_auroc(time_pred_combined_celltypei_sub$timerank, time_pred_combined_celltypei_sub$Cavebear)
          names(Cavebear_mat) <- c('gap', 'timerank', 'Cavebear')
          upperbaseline_mat <- calc_pairwise_auroc(time_pred_combined_celltypei_sub$timerank, time_pred_combined_celltypei_sub$upperbaseline)
          names(upperbaseline_mat) <- c('gap', 'timerank', 'UpperBaseline')
          monocle_mat <- calc_pairwise_auroc(time_pred_combined_celltypei_sub$timerank, time_pred_combined_celltypei_sub$monocle)
          names(monocle_mat) <- c('gap', 'timerank', 'Monocle')
          Cavebear_mat <- merge(Cavebear_mat, monocle_mat, by=c('gap', 'timerank'))
          Cavebear_mat <- merge(Cavebear_mat, upperbaseline_mat, by=c('gap', 'timerank'))
          Cavebear_mat$major_trajectory <- celltype_sub_i
          Cavebear_mat$celltype_sub_cate <- celltype_sub_cate
          auc_mat <- rbind(auc_mat, Cavebear_mat)
        }
      }
    }
  }

  if (FALSE){
    # only for presentation
    ## Cavebear only
    p_compare_meanline_bearonly  <- ggplot(auc_mat_melt_stat[auc_mat_melt_stat$method=='Cavebear',], aes(x=gap, y=precision, colour=method)) + facet_wrap(~major_trajectory) +
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
    
    ggsave(p_compare_meanline_bearonly, file=paste0(plot_dir, 'compare_pseudotime_alltime_', monocle_train_ver, 'bear_byct_barplot.png'), width=12, height=6)
  }
}

calc_pairwise_auroc_species <- function(time_pred_cavebear, cavebear_input_dir, plot_dir, data_ver){
  print(plot_dir)
  auc_mat <- c()
  if (data_ver=='human'){
    time_pred_cavebear$age_broad <- 1
    time_pred_cavebear[age>5]$age_broad <- 2
    time_pred_cavebear[age>5.5]$age_broad <- 3
    time_pred_cavebear$timerank <- dense_rank(time_pred_cavebear$age_broad)
  }else{
    time_pred_cavebear$timerank <- dense_rank(time_pred_cavebear$age)
  }
  
  ## calculate pair-wise AUROC on each cell type. try different cell type level to test the robustness of Cavebear
  for (celltype_sub_cate in c('predicted_cell_type')){
    time_pred_cavebear$celltype_sub_cate <- time_pred_cavebear[[celltype_sub_cate]]
    for (celltype_sub_i in unique(time_pred_cavebear$celltype_sub_cate)){
      time_pred_combined_celltypei_sub <- time_pred_cavebear[celltype_sub_cate==celltype_sub_i]
      if (length(unique(time_pred_combined_celltypei_sub$timerank))>2){
        Cavebear_mat <- calc_pairwise_auroc(time_pred_combined_celltypei_sub$timerank, time_pred_combined_celltypei_sub$Cavebear)
        names(Cavebear_mat) <- c('gap', 'timerank', 'Cavebear')
        Cavebear_mat$celltype_sub_cate <- celltype_sub_i
        auc_mat <- rbind(auc_mat, Cavebear_mat)
      }
    }
  }
  
  ## plot using mean and standard error
  auc_mat_melt <- melt(auc_mat, id.vars = c('gap', 'celltype_sub_cate'), measure.vars = c('Cavebear'), variable.name = c('method'), value.name = "precision")
  auc_mat_melt$gap <- factor(auc_mat_melt$gap, levels=unique(auc_mat_melt$gap)[order(unique(auc_mat_melt$gap))])
  auc_mat_melt_stat <- summarySE(auc_mat_melt, measurevar="precision", groupvars=c("gap", 'celltype_sub_cate', "method"))
  p_compare_meanline_bearonly  <- ggplot(auc_mat_melt_stat[auc_mat_melt_stat$method=='Cavebear',], aes(x=gap, y=precision, colour=method)) + facet_wrap(~celltype_sub_cate, ncol=5) +
    geom_errorbar(aes(ymin=precision-se, ymax=precision+se), width=.1) +
    geom_line() +
    geom_point()+
    xlab('Time gap') +
    ylab('Accuracy') +
    geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
    #theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
    theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
    theme(axis.text.x = element_text(angle = 75, vjust = 0.5))+ 
    guides(color = "none")
  print(unique(time_pred_cavebear$celltype_sub_cate))
  if (length(unique(time_pred_cavebear$celltype_sub_cate))>1){
    ggsave(p_compare_meanline_bearonly, file=paste0(plot_dir, 'compare_pseudotime_alltime_', data_ver, '_bear_byct_barplot.png'), width=12, height=6)
  }else{
    ggsave(p_compare_meanline_bearonly, file=paste0(plot_dir, 'compare_pseudotime_alltime_', data_ver, '_bear_byct_barplot.png'), width=4, height=2.8)
  }
  
}


## ==============================================================
## compare pairwise AUROC, trained and validated across all time points in target species (zebrafish)
## ==============================================================
if (target_species == 'zebrafish'){
  cavebear_input_dir <- '/proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-23/'
  monocle_input_dir <- '/proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-27monocle/'
  plot_dir <- '/proj/ranz_lab/users/nrittenhouse/proj/2025_nathouse_cs-pseudotime/figures/'
  best_icebear_para <- '_0.001_3_25' #best icebear parameters

  ## load Cavebear prediction
  prefix <- paste0(cavebear_input_dir, 'cvae_pytorch_best_model', best_icebear_para, '_')
  suffix <- '_mouse_pred_timezebrafish.txt'
  time_pred_cavebear <- load_existing_prediction(prefix, suffix)
  names(time_pred_cavebear)[ncol(time_pred_cavebear)] <- 'Cavebear'

  # focus on large cell trajectories (celltype_large), otherwise we may not have enough cells that cover many time points/return robust eval
  celltype_large <- names(table(time_pred_cavebear$major_trajectory))[table(time_pred_cavebear$major_trajectory)>10000] 

  ## load upper baseline
  prefix <- paste0(cavebear_input_dir, 'cvae_pytorch_best_model', best_icebear_para, '_')
  suffix <- '_zebrafish_pred_timezebrafish.txt'
  time_pred_upperbaseline <- load_existing_prediction(prefix, suffix)
  names(time_pred_upperbaseline)[ncol(time_pred_upperbaseline)] <- 'upperbaseline'
  time_pred_cavebear <- merge(time_pred_cavebear, time_pred_upperbaseline[, c('cell_unique', 'upperbaseline'), with=FALSE], by='cell_unique')

  # load monocle 2 prediction
  # because of memory constraint, monocle 2 is run separately on each cell trajectories. 

  ## for each cell trajectory, evaluate pseudotime prediction against actual time labels, compare between Cavebear, upperbaseline, and Monocle 2
  monocle_train_ver <- 'indiv_' #monocle is trained on each major cell trajectory separately (for memory's sake)

  compare_pairwise_auroc(celltype_large, time_pred_cavebear, monocle_input_dir, monocle_train_ver, cavebear_input_dir, plot_dir) # Figure 2B - old, failed to include two cell types & only compared with monocle 
  compare_pairwise_auroc(celltype_large, time_pred_cavebear, monocle_input_dir, monocle_train_ver, cavebear_input_dir, plot_dir) # Figure 2B - old, failed to include two cell types & only compared with monocle 
  #= Figure 2B: /proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-23/plot/compare_pseudotime_alltime_indiv_major_trajectoryupperbaseline_byct_barplot.pdf
  
  compare_pairwise_auroc_moremethods(celltype_large, time_pred_cavebear, monocle_input_dir, monocle_train_ver, cavebear_input_dir, palantir_input_dir, plot_dir) # Figure 2B - updated, all major cell trajectory and compared with multiple methods
  #= Figure 2B:/proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-23/plot/compare_pseudotime_alltime_indiv__all_methods_major_trajectory.pdf

  ## check small trajectories:
  celltype_small <- names(table(time_pred_cavebear$major_trajectory))[table(time_pred_cavebear$major_trajectory)<10000] 
  calc_pairwise_auroc_species(time_pred_cavebear[major_trajectory %in% celltype_small & !is.na(age)], cavebear_input_dir, plot_dir, 'zebrafish_smalltrajectory')
  #= decided to not show / show in supplementary: /proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-23/plot/compare_pseudotime_alltime_indiv_cell_type_subupperbaseline_byct_barplot.pdf
  
  ## ==============================================================
  ## compare correlation with upperbaseline, trained and validated across one time points in target species (zebrafish)
  ## ==============================================================
  cavebear_input_dir <- '/proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-27/'
  monocle_input_dir <- '/proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-27monocle/'
  plot_dir <- '/proj/ranz_lab/users/nrittenhouse/proj/2025_nathouse_cs-pseudotime/figures/'

  # start with upperbaseline
  #cell_annot_i <- time_pred_upperbaseline

  #based on VAE trained on mouse and one time point in zebrafish and predict on zebrafish
  timepoint_list <- c(18, 30, 42, 96)
  timepoint_list <- unique(time_pred_cavebear$age)
  #timepoint_list <- timepoint_list[!is.na(timepoint_list)]
  compare_singletime_correlation(celltype_large, time_pred_upperbaseline, monocle_input_dir, monocle_train_ver, best_icebear_para, timepoint_list, cavebear_input_dir, plot_dir)
  #= Figure 2D: /proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-27/plot/compare_pseudotime_pertime_indiv_abs_monocle2_byct.pdf
  #= Figure 2E: /proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-27/plot/compare_pseudotime_pertime_indiv_abs_monocle2_bytime.pdf
}