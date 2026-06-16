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

target_species = args[1] # options: zebrafish, mouse_EB, Cao, human, livercancer, cs_organoid_nosub

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
  
  ## plot using mean and standard error
  ## because Monocle is inferred on major trajectory level, we help Monocle by swapping time sign for each trajectory that is anti-correlated with real time labels.
  auc_mat[, monocle_sign := sign(mean(Monocle)-0.5), by=major_trajectory]
  auc_mat$Monocle_swap = auc_mat$Monocle
  auc_mat[monocle_sign < 0]$Monocle_swap = 1- auc_mat[monocle_sign < 0]$Monocle
  auc_mat_melt <- melt(auc_mat, id.vars = c('major_trajectory','gap', 'celltype_sub_cate'), measure.vars = c('Monocle', 'Monocle_swap', 'Cavebear', 'UpperBaseline'), variable.name = c('method'), value.name = "precision")
  auc_mat_melt$gap <- factor(auc_mat_melt$gap, levels=unique(auc_mat_melt$gap)[order(unique(auc_mat_melt$gap))])
  auc_mat_melt_stat <- summarySE(auc_mat_melt, measurevar="precision", groupvars=c("major_trajectory","gap", 'celltype_sub_cate', "method"))
  auc_mat_melt_stat$method <- factor(auc_mat_melt_stat$method, levels=c('Cavebear', 'Monocle', 'Monocle_swap', 'UpperBaseline'))

  
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
    
    ## compare Cavebear with original monocle 2 (sign not swapped)
    p_compare_meanline_ori  <- ggplot(auc_mat_melt_stat[auc_mat_melt_stat$method %in% c('Monocle', 'Cavebear'),], aes(x=gap, y=precision, colour=method)) + facet_wrap(~major_trajectory) +
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
    
    ggsave(p_compare_meanline_ori, file=paste0(plot_dir, 'compare_pseudotime_alltime_', monocle_train_ver, 'monocle2_byct_barplot.png'), width=12, height=6)
    
  }
  
  ## compare Cavebear with sign-swapped monocle 2
  for (celltype_sub_cate in c('major_trajectory', 'predicted_cell_type', 'cell_type_broad', 'cell_type_sub')){
    p_compare_meanline_swap  <- ggplot(auc_mat_melt_stat[auc_mat_melt_stat$method %in% c('Monocle_swap', 'Cavebear') & auc_mat_melt_stat$celltype_sub_cate==celltype_sub_cate,], aes(x=gap, y=precision, colour=method)) + facet_wrap(~major_trajectory, ncol=5) +
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
    
    ggsave(p_compare_meanline_swap, file=paste0(plot_dir, 'compare_pseudotime_alltime_', monocle_train_ver, '_', celltype_sub_cate, 'swap_monocle2_byct_barplot.png'), width=12, height=6)
    
    p_compare_meanline_upper  <- ggplot(auc_mat_melt_stat[auc_mat_melt_stat$method %in% c('UpperBaseline', 'Cavebear', 'Monocle_swap') & auc_mat_melt_stat$celltype_sub_cate==celltype_sub_cate,], aes(x=gap, y=precision, colour=method)) + facet_wrap(~major_trajectory) +
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
    
    ggsave(p_compare_meanline_upper, file=paste0(plot_dir, 'compare_pseudotime_alltime_', monocle_train_ver, celltype_sub_cate, 'upperbaseline_byct_barplot.png'), width=12, height=6)
    ggsave(p_compare_meanline_upper, file=paste0(plot_dir, 'compare_pseudotime_alltime_', monocle_train_ver, celltype_sub_cate, 'upperbaseline_byct_barplot.pdf'), width=12, height=6) # Figure 2B
  }
}


palantir_input_dir <- '/proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2026-05-20palantir/'
compare_pairwise_auroc_moremethods <- function(celltype_large, time_pred_cavebear, monocle_input_dir, monocle_train_ver, cavebear_input_dir, palantir_input_dir, plot_dir){
  print(plot_dir)
  auc_mat <- c()
  
  for (celltype_i in celltype_large){
    celltype_i_name <- gsub('/', '_', celltype_i)
    
    # ── load monocle ──────────────────────────────────────────────────────────
    monocle_file <- paste0(monocle_input_dir, 'monocle2/', 'pseudotime_', monocle_train_ver, celltype_i_name, '.txt')
    if (!file.exists(monocle_file)) next
    time_pred_monocle <- fread(monocle_file)
    
    # ── load palantir/dpt ─────────────────────────────────────────────────────
    palantir_file <- paste0(palantir_input_dir, 'pseudotime_indiv_', celltype_i_name, '.txt')
    has_palantir <- file.exists(palantir_file)
    if (has_palantir) time_pred_palantir <- fread(palantir_file)
    
    # ── merge cavebear + monocle ──────────────────────────────────────────────
    time_pred_combined <- merge(time_pred_cavebear, time_pred_monocle, by.x='cell_unique', by.y='cell')
    
    # ── merge palantir/dpt if available ──────────────────────────────────────
    if (has_palantir){
      time_pred_combined <- merge(
        time_pred_combined,
        time_pred_palantir[, .(cell_unique, palantir_pseudotime, dpt_pseudotime)],
        by='cell_unique', all.x=TRUE
      )
    } else {
      time_pred_combined$palantir_pseudotime <- NA_real_
      time_pred_combined$dpt_pseudotime      <- NA_real_
    }
    
    # ── filter to this cell type, remove bad rows ─────────────────────────────
    time_pred_combined_celltypei <- time_pred_combined[major_trajectory == celltype_i]
    time_pred_combined_celltypei <- time_pred_combined_celltypei[!is.infinite(monocle)]
    time_pred_combined_celltypei <- time_pred_combined_celltypei[!is.na(timepoint)]
    time_pred_combined_celltypei$timerank <- dense_rank(time_pred_combined_celltypei$age)
    
    # ── focus on major_trajectory level only ──────────────────────────────────
    for (celltype_sub_i in unique(time_pred_combined_celltypei$major_trajectory)){
      sub <- time_pred_combined_celltypei[major_trajectory == celltype_sub_i]
      
      Cavebear_mat      <- calc_pairwise_auroc(sub$timerank, sub$Cavebear)
      upperbaseline_mat <- calc_pairwise_auroc(sub$timerank, sub$upperbaseline)
      monocle_mat       <- calc_pairwise_auroc(sub$timerank, sub$monocle)
      palantir_mat      <- calc_pairwise_auroc(sub$timerank, sub$palantir_pseudotime)
      dpt_mat           <- calc_pairwise_auroc(sub$timerank, sub$dpt_pseudotime)
      
      names(Cavebear_mat)      <- c('gap', 'timerank', 'Cavebear')
      names(upperbaseline_mat) <- c('gap', 'timerank', 'UpperBaseline')
      names(monocle_mat)       <- c('gap', 'timerank', 'Monocle')
      names(palantir_mat)      <- c('gap', 'timerank', 'Palantir')
      names(dpt_mat)           <- c('gap', 'timerank', 'DPT')
      
      merged <- Reduce(function(a, b) merge(a, b, by=c('gap', 'timerank')),
                       list(Cavebear_mat, monocle_mat, palantir_mat, dpt_mat, upperbaseline_mat))
      merged$major_trajectory  <- celltype_sub_i
      merged$celltype_sub_cate <- 'major_trajectory'
      auc_mat <- rbind(auc_mat, merged)
    }
  }
  
  # ── sign-swap for all pseudotime methods ──────────────────────────────────
  auc_mat <- as.data.table(auc_mat)
  
  for (method_col in c('Monocle', 'Palantir', 'DPT')){
    swap_col <- paste0(method_col, '_swap')
    sign_col <- paste0(method_col, '_sign')
    auc_mat[, (sign_col) := sign(mean(get(method_col), na.rm=TRUE) - 0.5), by=major_trajectory]
    auc_mat[, (swap_col) := get(method_col)]
    auc_mat[get(sign_col) < 0, (swap_col) := 1 - get(method_col)]
  }
  # ── Rank major-trajectories by auc for Cavebear ──────────────────────────────────
  # rank by mean auc for all time gaps per a cell type
  rank_by_mean <- auc_mat[,
    .(mean_cavebear = mean(Cavebear)),
    by = major_trajectory
  ][order(-mean_cavebear)]
  
  # set factors to change order of major trajectories based on rank
  ordered_levels <- rank_by_mean[order(-mean_cavebear), major_trajectory]
  auc_mat[, major_trajectory := factor(major_trajectory, levels = ordered_levels)]

  # ── melt and plot ─────────────────────────────────────────────────────────
  measure_vars <- c('Cavebear', 'UpperBaseline',
                    'Monocle_swap', 'Palantir_swap', 'DPT_swap')

  auc_mat_melt <- melt(
    auc_mat,
    id.vars      = c('major_trajectory', 'gap', 'celltype_sub_cate'),
    measure.vars = measure_vars,
    variable.name = 'method',
    value.name    = 'precision'
  )
  auc_mat_melt$gap <- factor(auc_mat_melt$gap,
                             levels = unique(auc_mat_melt$gap)[order(unique(auc_mat_melt$gap))])
  
  auc_mat_melt_stat <- summarySE(auc_mat_melt,
                                 measurevar = 'precision',
                                 groupvars  = c('major_trajectory', 'gap', 'celltype_sub_cate', 'method'))
  
  method_levels <- c('Cavebear', 'Palantir_swap', 'DPT_swap', 'Monocle_swap', 'UpperBaseline')
  method_colors <- c('Cavebear'      = '#E41A1C',
                     'Palantir_swap' = '#377EB8',
                     'DPT_swap'      = '#4DAF4A',
                     'Monocle_swap'  = '#FF7F00',
                     'UpperBaseline' = '#984EA3')
  
  auc_mat_melt_stat$method <- factor(auc_mat_melt_stat$method, levels = method_levels)
  
  p_all <- ggplot(auc_mat_melt_stat,
                  aes(x=gap, y=precision, colour=method, group=method)) +
    facet_wrap(~major_trajectory, ncol = 5) +
    geom_errorbar(aes(ymin=precision-se, ymax=precision+se), width=.1) +
    geom_line() +
    geom_point() +
    scale_colour_manual(values=method_colors) +
    xlab('Time gap') +
    ylab('Pairwise AUROC') +
    geom_hline(yintercept=0.5, linetype='dashed', color='grey') +
    theme(
      panel.background = element_rect(fill='white', colour='white'),
      panel.border     = element_rect(colour='black', fill=NA, size=0.8),
      axis.text        = element_text(size=10, colour='black'),
      axis.title       = element_text(size=13, colour='black'),
      axis.text.x      = element_text(angle=75, vjust=0.5),
      legend.title     = element_blank()
    )
  
  ggsave(p_all,
         file  = paste0(plot_dir, 'compare_pseudotime_alltime_', monocle_train_ver, '_all_methods_major_trajectory.png'),
         width = 20, height = 12)
  ggsave(p_all,
         file  = paste0(plot_dir, 'compare_pseudotime_alltime_', monocle_train_ver, '_all_methods_major_trajectory.pdf'),
         width = 20, height = 12)
  
  invisible(auc_mat_melt_stat)
}

compare_pairwise_auroc_EB <- function(time_pred_cavebear, monocle_input_dir, monocle_train_ver,cavebear_input_dir, plot_dir){
  print(plot_dir)
  auc_mat <- c()

  filename <- paste0(monocle_input_dir, 'pseudotime_', monocle_train_ver, '.txt')
  if (file.exists(filename)) {
    time_pred_monocle <- fread(filename)
    time_pred_cavebear_monocle <- merge(time_pred_cavebear, time_pred_monocle, by = 'cell')
    time_pred_cavebear_monocle <- time_pred_cavebear_monocle[!is.infinite(monocle)]
    time_pred_cavebear_monocle$timerank <- dense_rank(time_pred_cavebear_monocle$mouse_age)
    print(unique(time_pred_cavebear_monocle$timerank))

    Cavebear_mat <- calc_pairwise_auroc(time_pred_cavebear_monocle$timerank, time_pred_cavebear_monocle$Cavebear)
    names(Cavebear_mat) <- c('gap', 'timerank', 'Cavebear')
    upperbaseline_mat <- calc_pairwise_auroc(time_pred_cavebear_monocle$timerank, time_pred_cavebear_monocle$upperbaseline)
    names(upperbaseline_mat) <- c('gap', 'timerank', 'UpperBaseline')
    monocle_mat <- calc_pairwise_auroc(time_pred_cavebear_monocle$timerank, time_pred_cavebear_monocle$monocle)
    names(monocle_mat) <- c('gap', 'timerank', 'Monocle')
    Cavebear_mat <- merge(Cavebear_mat, monocle_mat, by=c('gap', 'timerank'))
    Cavebear_mat <- merge(Cavebear_mat, upperbaseline_mat, by=c('gap', 'timerank'))
    auc_mat <- Cavebear_mat
  }

  ## plot using mean and standard error
  ## because Monocle is inferred on major trajectory level, we help Monocle by swapping time sign for each trajectory that is anti-correlated with real time labels.
  auc_mat[, monocle_sign := sign(mean(Monocle)-0.5)]
  auc_mat$Monocle_swap = auc_mat$Monocle
  auc_mat[monocle_sign < 0]$Monocle_swap = 1- auc_mat[monocle_sign < 0]$Monocle
  auc_mat_melt <- melt(auc_mat, id.vars = c('gap'), measure.vars = c('Monocle', 'Monocle_swap', 'Cavebear', 'UpperBaseline'), variable.name = c('method'), value.name = "precision")
  auc_mat_melt$gap <- factor(auc_mat_melt$gap, levels=unique(auc_mat_melt$gap)[order(unique(auc_mat_melt$gap))])
  auc_mat_melt_stat <- summarySE(auc_mat_melt, measurevar="precision", groupvars=c("gap", "method"))
  auc_mat_melt_stat$method <- factor(auc_mat_melt_stat$method, levels=c('Cavebear', 'Monocle', 'Monocle_swap', 'UpperBaseline'))
  
  if (TRUE){
    # only for presentation
    ## Cavebear only
    p_compare_meanline_bearonly  <- ggplot(auc_mat_melt_stat[auc_mat_melt_stat$method=='Cavebear',], aes(x=gap, y=precision, colour=method)) +
      geom_errorbar(aes(ymin=precision-se, ymax=precision+se), width=.1) +
      geom_line() +
      geom_point()+
      xlab('Time gap') +
      ylab('Accuracy') +
      geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
      #theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
      theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
      theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
      theme(axis.text.x = element_text(angle = 75, vjust = 0.5))
    
    ggsave(p_compare_meanline_bearonly, file=paste0(plot_dir, 'compare_pseudotime_alltime_', monocle_train_ver, 'bear_barplot.png'), width=12, height=6, create.dir = TRUE)
    
    ## compare Cavebear with original monocle 2 (sign not swapped)
    p_compare_meanline_ori  <- ggplot(auc_mat_melt_stat[auc_mat_melt_stat$method %in% c('Monocle', 'Cavebear'),], aes(x=gap, y=precision, colour=method)) +
      geom_errorbar(aes(ymin=precision-se, ymax=precision+se), width=.1) +
      geom_line() +
      geom_point()+
      xlab('Time gap') +
      ylab('Accuracy') +
      geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
      #theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
      theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
      theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
      theme(axis.text.x = element_text(angle = 75, vjust = 0.5))
    
    ggsave(p_compare_meanline_ori, file=paste0(plot_dir, 'compare_pseudotime_alltime_', monocle_train_ver, 'monocle2_barplot.png'), width=12, height=6, create.dir = TRUE)
    
  }
  
  ## compare Cavebear with sign-swapped monocle 2
  p_compare_meanline_swap  <- ggplot(auc_mat_melt_stat[auc_mat_melt_stat$method %in% c('Monocle_swap', 'Cavebear'),], aes(x=gap, y=precision, colour=method)) +
    geom_errorbar(aes(ymin=precision-se, ymax=precision+se), width=.1) +
    geom_line() +
    geom_point()+
    xlab('Time gap') +
    ylab('Accuracy') +
    geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
    #theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
    theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
    theme(axis.text.x = element_text(angle = 75, vjust = 0.5))
    
  ggsave(p_compare_meanline_swap, file=paste0(plot_dir, 'compare_pseudotime_alltime_', monocle_train_ver, 'swap_monocle2_barplot.png'), width=12, height=6,create.dir = TRUE)
    
  p_compare_meanline_upper  <- ggplot(auc_mat_melt_stat[auc_mat_melt_stat$method %in% c('UpperBaseline', 'Cavebear', 'Monocle_swap'),], aes(x=gap, y=precision, colour=method)) +
    geom_errorbar(aes(ymin=precision-se, ymax=precision+se), width=.1) +
    geom_line() +
    geom_point()+
    xlab('Time gap') +
    ylab('Accuracy') +
    geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
    #theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
    theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
    theme(axis.text.x = element_text(angle = 75, vjust = 0.5))
    
    ggsave(p_compare_meanline_upper, file=paste0(plot_dir, 'compare_pseudotime_alltime_', monocle_train_ver, 'upperbaseline_MonocleSwap_barplot.png'), width=12, height=6,create.dir = TRUE)

    p_compare_meanline_upper  <- ggplot(auc_mat_melt_stat[auc_mat_melt_stat$method %in% c('UpperBaseline', 'Cavebear', 'Monocle'),], aes(x=gap, y=precision, colour=method)) +
    geom_errorbar(aes(ymin=precision-se, ymax=precision+se), width=.1) +
    geom_line() +
    geom_point()+
    xlab('Time gap') +
    ylab('Accuracy') +
    geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
    #theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
    theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
    theme(axis.text.x = element_text(angle = 75, vjust = 0.5))
    
    ggsave(p_compare_meanline_upper, file=paste0(plot_dir, 'compare_pseudotime_alltime_', monocle_train_ver, 'upperbaseline_Monocle_barplot.png'), width=12, height=6)
} 


plot_prediction_distribution_violin <- function(data, plot_dir, time_label, method){
  levels <- unique(data[[time_label]])
  data <- data %>% rename("timepoint" = all_of(time_label))

  sample_size = data %>% group_by(timepoint) %>% summarize(num=n(), .groups = "drop")
  head(sample_size)

  # Plot
  data <- data %>% left_join(sample_size, by = "timepoint") %>% mutate(timepoint = factor(timepoint, levels=levels), myaxis = paste0(timepoint, "\n", "n=", num))
  violin_plot <- ggplot(data, aes(x=timepoint, y=.data[[method]], fill = timepoint)) +
    geom_violin(width=1.4) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    theme(legend.position="none", plot.title = element_text(size=11)) +
    xlab("")
  ggsave(violin_plot, file=paste0(plot_dir, 'compare_pseudotime_alltime_',method,'_violin.png'), width=4, height=2.8)
}


plot_prediction_distribution_violin_humansingleome <- function(data, plot_dir, time_label, method, model_name){
  levels <- unique(data[[time_label]])
  data <- data %>% rename("timepoint" = all_of(time_label))
  
  sample_size = data %>% group_by(timepoint) %>% summarize(num=n(), .groups = "drop")
  head(sample_size)
  
  # Plot
  data <- data %>% left_join(sample_size, by = "timepoint") %>% mutate(timepoint = factor(timepoint, levels=levels), myaxis = paste0(timepoint, "\n", "n=", num))
  data$timepoint <- factor(data$timepoint, levels=unique(data$timepoint)[order(as.numeric(as.character(unique(data$timepoint))))])
  violin_plot <- ggplot(data, aes(x=timepoint, y=.data[[method]])) +
    geom_violin(width=1.4) +
    geom_boxplot(width=0.1, notch=T) +
    theme(legend.position="none", plot.title = element_text(size=11)) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
          plot.background  = element_rect(fill = "white", colour = "white"),
          panel.grid       = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
    theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  
  ggsave(violin_plot, file=paste0(plot_dir, 'compare_pseudotime_alltime_', model_name, '_violinplot.png'), width=4, height=2.8)
  ggsave(violin_plot, file=paste0(plot_dir, 'compare_pseudotime_alltime_', model_name, '_violinplot.pdf'), width=4, height=2.8) # Figure 4A1
}


compare_singletime_correlation <- function(celltype_large, time_pred_upperbaseline, monocle_input_dir, monocle_train_ver, best_icebear_para, timepoint_list, cavebear_input_dir, plot_dir){
  cor_mat_combined <- c()
  for (timepoint_i in timepoint_list){
    print(timepoint_i)
    time_pred_upperbaseline_timepointi <- time_pred_upperbaseline[age==timepoint_i]
    # load monocle prediction
    monocle_filename <- paste0(monocle_input_dir, '/monocle2/pseudotime_', monocle_train_ver, timepoint_i, '.txt')
    match = 0
    if (file.exists(monocle_filename)){
      time_pred_monocle_timepointi <- fread(monocle_filename) #= monocle automatically removed some low quality cells. our prediction correspond to order of cell_annot_i_time, but need to merge with monocle after by cell id
      time_pred_monocle_timepointi$Monocle <- time_pred_monocle_timepointi$monocle
      # load cavebear prediction - with best time predictor
      prefix <- paste0(cavebear_input_dir, 'cvae_pytorch_best_model', timepoint_i, best_icebear_para, '_')
      suffix <- '_mouse_pred_timezebrafish.txt'
      time_pred_cavebear_timepointi <- load_existing_prediction(prefix, suffix)
      if (!is.null(time_pred_cavebear_timepointi)){ # if both monocle and cavebear prediction exists for that timepoint_i
        time_pred_cavebear_timepointi$Cavebear <- time_pred_cavebear_timepointi$pred_time
        time_pred_combined_timepointi <- merge(time_pred_upperbaseline_timepointi, time_pred_cavebear_timepointi[,c('cell_unique', 'Cavebear'), with=FALSE], by='cell_unique')
        time_pred_combined_timepointi <- merge(time_pred_combined_timepointi, time_pred_monocle_timepointi[,c('cell', 'Monocle'), with=FALSE], by.x='cell_unique', by.y = 'cell')
        cor_mat <- c()
        for (celltype_i in celltype_large){
          print(celltype_i)
          time_pred_combined_timepointi_celltypei <- time_pred_combined_timepointi[major_trajectory==celltype_i]
          time_pred_combined_timepointi_celltypei <- time_pred_combined_timepointi_celltypei[!is.infinite(Monocle)]
          scor_i <- cor(as.numeric(time_pred_combined_timepointi_celltypei$upperbaseline), time_pred_combined_timepointi_celltypei$Monocle, method='spearman')
          kcor_i <- kendall.tau(as.numeric(time_pred_combined_timepointi_celltypei$upperbaseline), time_pred_combined_timepointi_celltypei$Monocle, exact = FALSE, max.n = 3000)
          cor_mat <- rbind(cor_mat, c(celltype_i, scor_i, kcor_i, 'Monocle'))
          scor_i <- cor(as.numeric(time_pred_combined_timepointi_celltypei$upperbaseline), time_pred_combined_timepointi_celltypei$Cavebear, method='spearman')
          kcor_i <- kendall.tau(as.numeric(time_pred_combined_timepointi_celltypei$upperbaseline), time_pred_combined_timepointi_celltypei$Cavebear, exact = FALSE, max.n = 3000)
          cor_mat <- rbind(cor_mat, c(celltype_i, scor_i, kcor_i, 'Cavebear'))
        }
        cor_mat <- as.data.table(cor_mat)
        names(cor_mat) <- c('major_trajectory', 'Spearman', 'kendall.tau', 'method')
        cor_mat$time <- timepoint_i
        
        cor_mat_combined <- rbind(cor_mat_combined, cor_mat)
      }
      
    }
  }
  
  cor_mat_combined$Spearman <- as.numeric(cor_mat_combined$Spearman)
  cor_mat_combined_pairwise <- cor_mat_combined[method=='Cavebear']
  cor_mat_combined_pairwise$Cavebear <- cor_mat_combined_pairwise$Spearman
  cor_mat_combined_pairwise$Monocle <- cor_mat_combined[method!='Cavebear']$Spearman
  
  p_compare <- ggplot(cor_mat_combined_pairwise, aes(Monocle, Cavebear))+
    ggtitle(paste('Spearman correlation')) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    geom_point(size=2.5, stroke = 0, color='#3182bd')+ geom_abline(slope=1,intercept=0, color='#636363',linetype="dotted", size=1.2) +
    theme(legend.position = "none",axis.ticks = element_blank())+ theme(axis.text=element_text(size=12,colour="black"))+
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
    theme(axis.text=element_text(size=12,colour="black"), axis.title=element_text(size=22, colour='black'))
  #p_compare
  cor_mat_combined_pairwise$Monocle_swap <- abs(cor_mat_combined_pairwise$Monocle)
  pval <- wilcox.test(cor_mat_combined_pairwise$Cavebear, cor_mat_combined_pairwise$Monocle, alternative='greater', paired=TRUE)$p.value
  
  p_compare_abs <- ggplot(cor_mat_combined_pairwise, aes(x=Monocle_swap, y=Cavebear, color=major_trajectory))+
    ggtitle(paste('Spearman correlation')) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    geom_point(size=2.5, stroke = 0)+ geom_abline(slope=1,intercept=0, color='#636363',linetype="dotted", size=1.2) +
    #theme(legend.position = "none",axis.ticks = element_blank())+ 
    theme(axis.text=element_text(size=12,colour="black"))+
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
    theme(axis.text=element_text(size=12,colour="black"), axis.title=element_text(size=22, colour='black')) +
    annotate(geom="text", x = max(cor_mat_combined_pairwise$Monocle_swap, na.rm = TRUE)-0.15, y = 0.1, label=paste0("P = ",formatC(pval, format = "e", digits = 2)),color="black", size=5)
  #p_compare_abs
  p_compare_abs_coltime <- ggplot(cor_mat_combined_pairwise, aes(x=Monocle_swap, y=Cavebear, color=time))+
    ggtitle(paste('Spearman correlation')) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    geom_point(size=2.5, stroke = 0)+ geom_abline(slope=1,intercept=0, color='#636363',linetype="dotted", size=1.2) +
    #theme(legend.position = "none",axis.ticks = element_blank())+ 
    theme(axis.text=element_text(size=12,colour="black"))+
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
    theme(axis.text=element_text(size=12,colour="black"), axis.title=element_text(size=22, colour='black'))
  #p_compare_abs_coltime
  
  ggsave(p_compare, file=paste0(plot_dir, 'compare_pseudotime_pertime_', monocle_train_ver, 'monocle2.png'), width=4.1, height=4.5)
  ggsave(p_compare_abs, file=paste0(plot_dir, 'compare_pseudotime_pertime_', monocle_train_ver, 'abs_monocle2_byct.png'), width=8.1, height=4.5)
  ggsave(p_compare_abs, file=paste0(plot_dir, 'compare_pseudotime_pertime_', monocle_train_ver, 'abs_monocle2_byct.pdf'), width=8.1, height=4.5) # Figure 2C
  ggsave(p_compare_abs_coltime, file=paste0(plot_dir, 'compare_pseudotime_pertime_', monocle_train_ver, 'abs_monocle2_bytime.png'), width=5, height=4.5)
  ggsave(p_compare_abs_coltime, file=paste0(plot_dir, 'compare_pseudotime_pertime_', monocle_train_ver, 'abs_monocle2_bytime.pdf'), width=5, height=4.5) # Figure 2D
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
  ## because Monocle is inferred on major trajectory level, we help Monocle by swapping time sign for each trajectory that is anti-correlated with real time labels.
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
  #= Figure 2C: /proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-27/plot/compare_pseudotime_pertime_indiv_abs_monocle2_byct.pdf
  #= Figure 2D: /proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-27/plot/compare_pseudotime_pertime_indiv_abs_monocle2_bytime.pdf

  ## ==============================================================
  ## use mouse to predict human age (Cao et al.)
  ## ==============================================================
  ## load Cavebear prediction

  cavebear_input_dir <- '/proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-27/'
  plot_dir <- '/proj/ranz_lab/users/nrittenhouse/proj/2025_nathouse_cs-pseudotime/figures/'
  best_icebear_para <- '_0.001_3_25human' #best icebear parameters

  prefix <- paste0(cavebear_input_dir, 'cvae_pytorch_best_model', best_icebear_para, '_')
  suffix <- '_mouse_pred_timehuman.txt'
  time_pred_cavebear <- load_existing_prediction(prefix, suffix)
  names(time_pred_cavebear)[ncol(time_pred_cavebear)] <- 'Cavebear'
  #calc_pairwise_auroc_human(time_pred_cavebear, cavebear_input_dir, plot_dir, 'human')
  #= the result is not very convincing - also differ from previous prediction - need to optimize model, not sure if it's a good direction for this paper, can be in the multi-species label project
}


## ==============================================================
## compare pairwise AUROC, trained and validated across all time points in target species (mouse_EB) for in vitro predictions
## ==============================================================
if (target_species == 'mouse_EB'){
  cavebear_input_dir <- '/proj/ranz_lab/users/nrittenhouse/proj/2025_nathouse_cs-pseudotime/results/2026-01-21/'
  monocle_input_dir <- '/proj/ranz_lab/users/nrittenhouse/proj/2025_nathouse_cs-pseudotime/results/monocle2_invitro_mouseEB/'
  plot_dir <- '/proj/ranz_lab/users/nrittenhouse/proj/2025_nathouse_cs-pseudotime/figures/'
  best_icebear_para <- '__0.001_3_25_Noble_EB_earlyEmbryonic_updatedBatch' #best icebear parameters
  time_label_bear <- 'timepoint'
  time_label_mon <- 'collection_time'

  ## load Cavebear prediction
  prefix <- paste0(cavebear_input_dir, 'cvae_pytorch_best_model', best_icebear_para, '_')
  suffix <- '_mouse_pred_time_mouse_EB.txt'
  time_pred_cavebear <- load_existing_prediction(prefix, suffix)
  names(time_pred_cavebear)[ncol(time_pred_cavebear)] <- 'Cavebear'

  ## load upper baseline
  prefix <- paste0(cavebear_input_dir, 'cvae_pytorch_best_model', best_icebear_para, '_')
  suffix <- '_mouse_EB_pred_time_mouse_EB.txt'
  time_pred_upperbaseline <- load_existing_prediction(prefix, suffix)
  names(time_pred_upperbaseline)[ncol(time_pred_upperbaseline)] <- 'upperbaseline'
  time_pred_cavebear <- merge(time_pred_cavebear, time_pred_upperbaseline[, c('cell', 'upperbaseline'), with=FALSE], by='cell')

  ## load monocle 2 prediction
  monocle_train_ver <- 'indiv_invitro_mouseEB'
  filename <- paste0(monocle_input_dir, 'pseudotime_', monocle_train_ver, '.txt')
  time_pred_monocle <- fread(filename)
  names(time_pred_monocle)[ncol(time_pred_monocle)] <- 'Monocle'

  #compare_pairwise_auroc_EB(time_pred_cavebear, monocle_input_dir, monocle_train_ver, cavebear_input_dir, plot_dir) 

  ## Create Violin Plot of predicted time for time points
  plot_prediction_distribution_violin(time_pred_cavebear, plot_dir, time_label_bear, "Cavebear")
  plot_prediction_distribution_violin(time_pred_monocle, plot_dir, time_label_mon, "Monocle")

}


## ==============================================================
## use mouse to predict human age (Cao et al.)
## ==============================================================

if (target_species == 'Cao'){
  cavebear_input_dir <- '/proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/results/2025-12-27/'
  plot_dir <- '/proj/ranz_lab/users/nrittenhouse/proj/2025_nathouse_cs-pseudotime/figures/'
  best_icebear_para <- '_0.001_3_25human' #best icebear parameters
  
  prefix <- paste0(cavebear_input_dir, 'cvae_pytorch_best_model', best_icebear_para, '_')
  suffix <- '_mouse_pred_timehuman.txt'
  time_pred_cavebear <- load_existing_prediction(prefix, suffix)
  names(time_pred_cavebear)[ncol(time_pred_cavebear)] <- 'Cavebear'
  calc_pairwise_auroc_human(time_pred_cavebear, cavebear_input_dir, plot_dir, 'human')
  #= the result is not very convincing - also differ from previous prediction - need to optimize model, not sure if it's a good direction for this paper, can be in the multi-species label project
}

## ==============================================================
## evaluate human singleome scRNA-seq cell time prediction per cell type using pairwise AUROC
## ==============================================================
if (target_species == 'human'){
  for (model_name in c('humanSingleome_mouseNeuro_dis50.0', 'humanSingleome_mouseNeuro_dis5.0', 'humanSingleome_mouseNeuro_dis10.0')){
    plot_dir <- '/proj/ranz_lab/users/nrittenhouse/proj/2025_nathouse_cs-pseudotime/figures/'
    time_label_bear <- 'age'
    
    ## load Cavebear prediction
    if (model_name == 'humanSingleome_mouseNeuro_dis50.0'){
      time_pred_cavebear <- fread(paste0('/proj/ranz_lab/users/nrittenhouse/proj/2025_nathouse_cs-pseudotime/results/2026-04-13/humanSingleome_mouseNeuro/0.001_3_25_batch/cvae_pytorch_disc_best_model__0.001_3_25_', model_name, '_human_pred_time_human.txt'))
    }else{
      time_pred_cavebear <- fread(paste0('/proj/ranz_lab/users/nrittenhouse/proj/2025_nathouse_cs-pseudotime/results/2026-04-10/humanSingleome_mouseNeuro/0.001_3_25_batch/cvae_pytorch_disc_best_model__0.001_3_25_', model_name, '_human_pred_time_human.txt'))
    }
    
    # optional: because pcw20 and pcw21 is too close, check if performance improves when we group them
    #time_pred_cavebear[sampleAge=='pcw20']$sampleAge <- 'pcw21'
    
    time_pred_cavebear[, age := as.numeric(sub("^pcw", "", sampleAge))]
    
    ## pairwise AUROC
    time_pred_cavebear$timerank <- dense_rank(time_pred_cavebear$age)
    time_pred_cavebear$Cavebear <- as.numeric(time_pred_cavebear$pred_time)
    
    auc_mat <- c()
    ## calculate pair-wise AUROC on each cell type. try different cell type level to test the robustness of Cavebear
    for (celltype_sub_cate in c('CellType_short')){
      time_pred_cavebear$celltype_sub_cate <- time_pred_cavebear[[celltype_sub_cate]]
      for (celltype_sub_i in unique(time_pred_cavebear$celltype_sub_cate)){
        time_pred_combined_celltypei_sub <- time_pred_cavebear[celltype_sub_cate==celltype_sub_i]
        # filter out ones with less than 100 cells
        time_pred_combined_celltypei_sub <- time_pred_combined_celltypei_sub[, if (.N >= 25) .SD, by = sampleAge]
        if (length(unique(time_pred_combined_celltypei_sub$age)) > 1 ) {
          print(celltype_sub_i)
          Cavebear_mat <- calc_pairwise_auroc(time_pred_combined_celltypei_sub$timerank, time_pred_combined_celltypei_sub$Cavebear)
          names(Cavebear_mat) <- c('gap', 'timerank', 'Cavebear')
          #Cavebear_mat$major_trajectory <- celltype_sub_i
          Cavebear_mat$celltype_sub_cate <- celltype_sub_i
          Cavebear_mat$ncell <- nrow(time_pred_combined_celltypei_sub)
          auc_mat <- rbind(auc_mat, Cavebear_mat)
        }
      }
    }
    auc_mat <- auc_mat[!is.na(Cavebear)]
    
    auc_mat_melt_stat <- summarySE(auc_mat, measurevar="Cavebear", groupvars=c("gap", "celltype_sub_cate"))
    p_compare_meanline_ori  <- ggplot(auc_mat_melt_stat, aes(x=gap, y=Cavebear)) + facet_wrap(~celltype_sub_cate, nrow=2) +
      geom_errorbar(aes(ymin=Cavebear-se, ymax=Cavebear+se), width=.1) +
      #geom_line() +
      geom_point()+
      xlab('Time gap') +
      ylab('Accuracy') +
      geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
      #theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
      theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
      theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
      scale_x_continuous(
        breaks = scales::breaks_width(1),
        labels = scales::label_number(accuracy = 1)
      )
    
    ggsave(p_compare_meanline_ori, file=paste0(plot_dir, 'compare_pseudotime_alltime_', model_name, '_barplot.png'), width=5.5, height=4, create.dir = TRUE)
    ggsave(p_compare_meanline_ori, file=paste0(plot_dir, 'compare_pseudotime_alltime_', model_name, '_barplot.pdf'), width=5.5, height=4, create.dir = TRUE) # Figure 4A2
    #= Figure 4A2: /proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/src/plot/compare_pseudotime_alltime_humanSingleome_mouseNeuro_dis10.0_barplot.pdf
    # TODO: please rerun with the best hyperparameter
    
    # overall violin plot
    plot_prediction_distribution_violin_humansingleome(time_pred_cavebear, plot_dir, time_label_bear, "Cavebear", model_name) # Figure 4A1
    #= Figure 4A1: /proj/ranz_lab/users/ranzhang/proj/2025_nathouse_cs-pseudotime/src/plot/compare_pseudotime_alltime_humanSingleome_mouseNeuro_dis10.0_barplot.pdf
    # TODO: please rerun with the best hyperparameter
  }
}


## ==============================================================
## evaluate human liver cancer time prediction
## ==============================================================
if (target_species == 'livercancer'){
  time_ver <- '2026-04-28liver'
  #time_ver <- '2026-05-18liverrescale'
  ## 1. besides the standard umap output, check umap of aligned liver fetal vs cancer samples:
  sim_url <- paste0('/proj/ranz_lab/users/ranzhang/proj/2025_ranzhang_cs-pseudotime/results/', time_ver, '/results/umaps/cvae_pytorch_best_model__0.001_3_25_cs_liver_updated_human')
  #sim_url <- paste0('/proj/ranz_lab/users/ranzhang/proj/2025_ranzhang_cs-pseudotime/results/', time_ver, '/results/umaps/cvae_pytorch_best_model__0.001_3_25_cs_liver_updated_dis4.0_human')
  #if (time_ver=='2026-04-28liver'){
  if (FALSE) {
    sim_url <- paste0('/proj/ranz_lab/users/ranzhang/proj/2025_ranzhang_cs-pseudotime/results/', time_ver, '/results/umaps/cvae_pytorch_best_model__0.0001_3_25_cs_liver_updated_human')
    embedding_mat_ori <- fread(paste0(sim_url, '_umap.txt'), sep='\t', header=F)
    embedding_mat <- fread(paste0(sim_url, '_umap.txt'), sep='\t', header=F)
    embedding_mat$V4 <- embedding_mat_ori$V4
    #embedding_mat[, (names(embedding_mat)) := lapply(.SD, function(x) gsub("^b'|'", "", x))]
    names(embedding_mat) <- c('UMAP_1', 'UMAP_2','batch', 'group', 'species')
    embedding_mat <- as.data.frame(embedding_mat)
    #embedding_mat$dataset <- factor(embedding_mat$dataset)
    embedding_mat$UMAP_1 <- as.numeric(embedding_mat$UMAP_1)
    embedding_mat$UMAP_2 <- as.numeric(embedding_mat$UMAP_2)
    embedding_mat$group <- factor(embedding_mat$group)
    embedding_mat$batch <- factor(embedding_mat$batch)
    embedding_mat$species <- factor(embedding_mat$species)
    embedding_mat <- as.data.table(embedding_mat)
    print(levels(embedding_mat$species))
    
    ## plot all, separation by batch, time or cell type
    dsize = 0.5
    transparency = 0.2
    resolution = 110
    
    p1_pred_combined_umap_ct_human <- ggplot(embedding_mat[species=='human'], aes(UMAP_1, UMAP_2, color = batch)) + facet_wrap(~group) +
      theme_classic() + theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
                              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
      geom_point(size=dsize, alpha = transparency) 
    
    p1_pred_combined_umap_ct_mouse <- ggplot(embedding_mat[species=='mouse' | batch=='liver'], aes(UMAP_1, UMAP_2, color = batch)) + facet_wrap(~group) +
      theme_classic() + theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
                              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
      geom_point(size=dsize, alpha = transparency) 
    
    p1_pred_combined_umap_ct_all <- ggplot(embedding_mat, aes(UMAP_1, UMAP_2, color = batch)) + facet_wrap(~group) +
      theme_classic() + theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
                              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
      geom_point(size=dsize, alpha = transparency) 
    
    unique(embedding_mat[species=='mouse']$group)
    
    png(paste0(sim_url, '_umap_group_colored_batch.png'), width = 2000, height = 1000, res=resolution)
    print(p1_pred_combined_umap_ct_all)
    dev.off()
    
    p1_pred_combined_umap_ct_fetal <- ggplot(embedding_mat[batch %in% c('fetal','fetal2'),], aes(UMAP_1, UMAP_2, color = species)) + facet_wrap(~group) +
      theme_classic() + theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
                              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
      geom_point(size=dsize, alpha = transparency) 
    
    p1_pred_combined_umap_ct_fetal
  }
  
  ## 2. check time prediction and if they differ between cancer vs normal
  plot_dir <- '/proj/ranz_lab/users/nrittenhouse/proj/2025_nathouse_cs-pseudotime/figures/'
  dir.create(plot_dir)
  time_label_bear <- 'age'
  
  ## load Cavebear prediction
  cell_annot <- fread('/proj/ranz_lab/users/ranzhang/proj/2025_ranzhang_cs-pseudotime/data/Cell2020Liver/cell_meta.csv')
  sample_annot <- fread('/proj/ranz_lab/users/ranzhang/proj/2025_ranzhang_cs-pseudotime/data/Cell2020Liver/sample_annot.txt')
  #method <- '0.0001_3_25_cs_liver_updated'
  method <- '0.001_3_25_cs_liver_updated'
  #method <- '0.001_3_25_cs_liver_updated_dis4.0'
  time_pred_cavebear <- fread(paste0('/proj/ranz_lab/users/ranzhang/proj/2025_ranzhang_cs-pseudotime/results/', time_ver, '/results/cvae_pytorch_best_model__', method, '_human_pred_time_human.txt'))
  # merge with fine-scaled cell type/groups
  time_pred_cavebear <- time_pred_cavebear[sampleID != 'HN']
  time_pred_cavebear <- time_pred_cavebear[batch=='liver']
  time_pred_cavebear$NTF <- factor(time_pred_cavebear$NTF, levels=unique(as.character(time_pred_cavebear$NTF)))
  time_pred_cavebear <- merge(time_pred_cavebear, cell_annot[, c('V1', 'CT_mye'), with=FALSE], by.x='cell_id', by.y='V1')
  
  trajectory_name <- 'major_trajectory' 
  trajectory_name <- 'CT_mye' # a more detailed cell annotation for myeloid
  time_pred_cavebear$trajectory <- time_pred_cavebear[[trajectory_name]]
  
  # filter out trajectory with very low number of cells
  if (trajectory_name=='CT_mye'){
    nmin_cell <- 5
  }else{
    nmin_cell <- 5
  }
  time_pred_cavebear <- time_pred_cavebear[, if (.N >= nmin_cell) .SD, by = .(sampleID, trajectory, NTF)]
  # step 1: retain CT_mye × sampleID combos that have both Tumor and Adj Normal
  time_pred_cavebear <- time_pred_cavebear[
    , if (all(c("Tumor", "Adj Normal") %in% NTF)) .SD,
    by = .(CT_mye, sampleID)
  ]
  
  # step 2: retain CT_mye with more than 5 unique sampleIDs
  time_pred_cavebear <- time_pred_cavebear[
    , if (uniqueN(sampleID) >= 3) .SD,
    by = CT_mye
  ]
  
  # verify
  time_pred_cavebear[, .(n_samples = uniqueN(sampleID), 
                         ntf_classes = paste(sort(unique(NTF)), collapse=',')),
                     by = CT_mye]
  
  ## load sample annotations to see if adjacent normal's prediction agrees with age:
  time_pred_cavebear <- merge(time_pred_cavebear, sample_annot[, c('Individual ID', 'HBV_Status (Clinical reports)', 'Gender', 'Age (2019)'), with=FALSE], by.x='sampleID', by.y='Individual ID')
  
  time_pred_cavebear$SampleAge <- factor(time_pred_cavebear$`Age (2019)`, 
                                         levels = unique(time_pred_cavebear$`Age (2019)`)[order(unique(time_pred_cavebear$`Age (2019)`))])
  
  # violin plot
  violin_plot_bysample <- ggplot(time_pred_cavebear, aes(x=sampleID, y=pred_time, fill = NTF)) + facet_wrap(~trajectory) +
    geom_violin(width=1.2, position = position_dodge(width = 0.8)) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
    theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  
  ggsave(violin_plot_bysample, file=paste0(plot_dir, 'compare_pseudotime_alltime_',method, trajectory_name, '_bysample_violin.png'), width=10, height=5)
  ggsave(violin_plot_bysample, file=paste0(plot_dir, 'compare_pseudotime_alltime_',method, trajectory_name, '_bysample_violin.pdf'), width=10, height=5) # Figure 4C
  # Figure 4C: /proj/ranz_lab/users/ranzhang/proj/2025_ranzhang_cs-pseudotime/results/2026-04-28liver/plot/compare_pseudotime_alltime_0.0001_3_25_cs_liver_updated_bysample_violin.pdf
  
  if (FALSE){
    ## check if there's difference by age
    violin_plot_byage <- ggplot(time_pred_cavebear, aes(x=SampleAge, y=pred_time, fill = NTF)) + facet_wrap(~trajectory) +
      geom_violin(width=1.2, position = position_dodge(width = 0.8)) +
      theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
      theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
    ggsave(violin_plot_byage, file=paste0(plot_dir, 'compare_pseudotime_alltime_',method,'_byage_violin.png'), width=10, height=2.8)
    
    ## check if there's sex difference
    violin_plot_byage <- ggplot(time_pred_cavebear[NTF=='Adj Normal'], aes(x=SampleAge, y=pred_time, fill = Gender)) + facet_wrap(~trajectory) +
      geom_violin(width=1.2, position = position_dodge(width = 0.8)) +
      theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
      theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
    
    ## check if there's influence from HBV_Status (Clinical reports)
    violin_plot_byage <- ggplot(time_pred_cavebear[NTF=='Adj Normal'], aes(x=SampleAge, y=pred_time, fill = `HBV_Status (Clinical reports)`)) + facet_wrap(~trajectory) +
      geom_violin(width=1.2, position = position_dodge(width = 0.8)) +
      theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
      theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
    
    ## check sub cell types
    violin_plot_bysample <- ggplot(time_pred_cavebear, aes(x=sampleID, y=pred_time, fill = NTF)) + facet_wrap(~trajectory) +
      geom_violin(width=1.2, position = position_dodge(width = 0.8)) +
      theme(panel.background = element_rect(fill = 'white', colour = 'white'), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8)) +
      theme(axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=13, colour='black')) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  }
  
  ## statistical test on whether endothelial and macrophage cells are shifting towards early development
  # take the mean per group, then do paired wilcoxon
  
  # mean pred_time per sampleID × group × NTF
  dt_mean <- time_pred_cavebear[
    , .(mean_pred_time = mean(pred_time, na.rm = TRUE)),
    by = .(sampleID, trajectory, NTF)
  ]
  
  # paired Wilcoxon test within each trajectory
  case_col <- "Tumor"
  control_col <- "Adj Normal"
  
  res <- dt_mean[
    , {
      wide <- dcast(.SD, sampleID ~ NTF, value.var = "mean_pred_time")
      
      if (all(c(case_col, control_col) %in% names(wide))) {
        
        cols <- c(case_col, control_col)
        complete <- wide[complete.cases(wide[, cols, with = FALSE])]
        
        if (nrow(complete) >= 2) {
          test <- wilcox.test(
            complete[[case_col]],
            complete[[control_col]],
            paired = TRUE,
            alternative = "less"
          )
          
          .(
            N_pairs = nrow(complete),
            mean_case = mean(complete[[case_col]], na.rm = TRUE),
            mean_control = mean(complete[[control_col]], na.rm = TRUE),
            mean_diff_case_minus_control =
              mean(complete[[case_col]] - complete[[control_col]], na.rm = TRUE),
            p_value = test$p.value
          )
        } else {
          .(
            N_pairs = nrow(complete),
            mean_case = NA_real_,
            mean_control = NA_real_,
            mean_diff_case_minus_control = NA_real_,
            p_value = NA_real_
          )
        }
        
      } else {
        .(
          N_pairs = NA_integer_,
          mean_case = NA_real_,
          mean_control = NA_real_,
          mean_diff_case_minus_control = NA_real_,
          p_value = NA_real_
        )
      }
    },
    by = trajectory
  ]
  
  res$FDR <- p.adjust(res$p_value, method = 'BH')
  View(res[order(FDR)])
  write.table(res, paste0(plot_dir, 'compare_pseudotime_alltime_',method, trajectory_name, '_stat.txt'), quote=F, row.names=F, col.names=T, sep='\t')
  ## let's check for cell cycle and stress status and correct for that in our prediction
}


## ==============================================================
## evaluate organoid pseudotime prediction (cs_organoid_nosub)
## ==============================================================
if (target_species == 'cs_organoid_nosub'){
  result_dir  <- '/proj/ranz_lab/users/ranzhang/proj/2025_ranzhang_cs-pseudotime/results/2026-04-05organoid_sub'
  data_dir <- '/proj/ranz_lab/users/ranzhang/proj/2025_ranzhang_cs-pseudotime/results/2026-04-01organoid'
  plot_dir <- '/proj/ranz_lab/users/nrittenhouse/proj/2025_nathouse_cs-pseudotime/figures/'
  dir.create(plot_dir, showWarnings = FALSE)
  
  dodge_narrow <- position_dodge(width = 0.8)
  theme_panel  <- theme(
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.text        = element_text(size = 10, colour = "black"),
    axis.title       = element_text(size = 13, colour = 'black'),
    axis.text.x      = element_text(angle = 45, vjust = 0.5)
  )
  
  # ── load data ──────────────────────────────────────────────────────────────
  pred       <- fread(file.path(result_dir, 'results/cvae_pytorch_best_model__0.001_3_25_cs_organoid_nosub_dis20.0_human_pred_time_human.txt'))
  cell_annot <- fread(file.path(data_dir, 'cs_organoid_nosub.csv'))
  
  pred[, day := sapply(strsplit(sampleNames, "_"), tail, 1)]
  
  pred[, DonorID := as.character(DonorID)]
  
  # ── overall pseudotime violin by day ─────────────────────────────────────
  p_pred_overall <- ggplot(pred[species=='human'], aes(y = pred_time, x = day)) +
    geom_violin(position = dodge_narrow, scale = "width") +
    geom_boxplot(width = 0.1, notch = TRUE, position = dodge_narrow,
                 outlier.shape = NA, color = "black") +
    ylab('predicted pseudotime') + xlab('') +
    theme_classic() + theme_panel
  ggsave(p_pred_overall, file = paste0(plot_dir, 'cs_organoid_nosub_pred_overall.png'), width = 6, height = 4)
  ggsave(p_pred_overall, file = paste0(plot_dir, 'cs_organoid_nosub_pred_overall.pdf'), width = 6, height = 4)
  
  # ── pseudotime by cell type ───────────────────────────────────────────────
  p_pred_by_celltype <- ggplot(pred[species=='human'], aes(y = pred_time, x = day)) +
    facet_wrap(~CellClass, ncol=4) +
    geom_violin(position = dodge_narrow, scale = "width") +
    geom_boxplot(width = 0.1, notch = TRUE, position = dodge_narrow,
                 outlier.shape = NA, color = "black") +
    ylab('predicted pseudotime') + xlab('') +
    theme_classic() + theme_panel +
    theme(legend.position = "bottom")
  ggsave(p_pred_by_celltype, file = paste0(plot_dir, 'cs_organoid_nosub_pred_celltype.png'), width = 20, height = 15, units = "cm", dpi = 600)
  ggsave(p_pred_by_celltype, file = paste0(plot_dir, 'cs_organoid_nosub_pred_celltype.pdf'), width = 20, height = 15, units = "cm") # Figure 3C
  # Figure 3C: /proj/ranz_lab/users/ranzhang/proj/2025_ranzhang_cs-pseudotime/results/2026-04-05organoid_sub/plot/cs_organoid_nosub_pred_celltype.pdf
  
  if (FALSE){
    # ── load individual-level phenotype ───────────────────────────────────────
    indiv_annot <- read_excel(paste0(result_dir, "/iviv_participant_info.xlsx"), sheet = 1)
    indiv_annot <- data.table(indiv_annot)
    indiv_annot[, DCCID        := as.character(DCCID)]
    indiv_annot[, ASD_Ever_Dx  := `ASD Ever Dx`]
    pred[indiv_annot, Sex         := i.Sex,         on = .(DonorID = DCCID)]
    pred[indiv_annot, Race        := i.Race,         on = .(DonorID = DCCID)]
    pred[indiv_annot, ASD_Ever_Dx := i.ASD_Ever_Dx, on = .(DonorID = DCCID)]
    
    # ── sex difference (ASD- only) ────────────────────────────────────────────
    p_pred_sex_group <- ggplot(pred[species=='human' & ASD_Ever_Dx=='ASD-'], aes(y = pred_time, x = day, fill = Sex, color = Sex)) +
      facet_wrap(~CellClass, ncol=4) +
      geom_violin(position = dodge_narrow, scale = "width", alpha = 0.6) +
      geom_boxplot(width = 0.1, notch = TRUE, position = dodge_narrow,
                   outlier.shape = NA, color = "black") +
      scale_fill_manual(values  = c("#fec44f", "#3182bd")) +
      scale_color_manual(values = c("#fec44f", "#3182bd")) +
      ylab('predicted pseudotime') + xlab('') +
      theme_classic() + theme_panel +
      theme(legend.position = "bottom")
    ggsave(p_pred_sex_group, file = paste0(plot_dir, 'cs_organoid_nosub_pred_sex_celltype.png'), width = 20, height = 15, units = "cm", dpi = 600)
    ggsave(p_pred_sex_group, file = paste0(plot_dir, 'cs_organoid_nosub_pred_sex_celltype.pdf'), width = 20, height = 15, units = "cm")
    
    # ── ASD difference (Male only) ────────────────────────────────────────────
    p_pred_asd_group <- ggplot(pred[species=='human' & Sex=='Male'], aes(y = pred_time, x = day, fill = ASD_Ever_Dx, color = ASD_Ever_Dx)) +
      facet_wrap(~CellClass) +
      geom_violin(position = dodge_narrow, scale = "width", alpha = 0.6) +
      geom_boxplot(width = 0.1, notch = TRUE, position = dodge_narrow,
                   outlier.shape = NA, color = "black") +
      scale_fill_manual(values  = c("#fec44f", "#3182bd")) +
      scale_color_manual(values = c("#fec44f", "#3182bd")) +
      ylab('predicted pseudotime') + xlab('') +
      theme_classic() + theme_panel +
      theme(legend.position = "bottom")
    ggsave(p_pred_asd_group, file = paste0(plot_dir, 'cs_organoid_nosub_pred_asd_celltype.png'), width = 20, height = 15, units = "cm", dpi = 600)
    ggsave(p_pred_asd_group, file = paste0(plot_dir, 'cs_organoid_nosub_pred_asd_celltype.pdf'), width = 20, height = 15, units = "cm")
    
    # ── statistical tests ─────────────────────────────────────────────────────
    # sex difference
    pred_forstats <- pred[species=='human' & ASD_Ever_Dx=='ASD-']
    pred_means    <- pred_forstats[, .(mean_pred_time = median(pred_time)), by = .(DonorID, day, Sex, CellClass, batch, batchID)]
    pred_means    <- pred_means[,   .(mean_pred_time = median(mean_pred_time)), by = .(DonorID, day, Sex, CellClass, batchID)]
    pval_mat <- c()
    for (day_i in unique(pred_means$day)){
      for (ct_i in unique(pred_means$CellClass)){
        pval <- tryCatch(
          wilcox.test(mean_pred_time ~ Sex, data = pred_means[day==day_i & CellClass==ct_i])$p.value,
          error = function(e) NA_real_)
        pval_mat <- rbind(pval_mat, c(ct_i, day_i, pval))
      }
    }
    pval_mat_sex        <- as.data.table(pval_mat)
    names(pval_mat_sex) <- c('CellClass', 'day', 'pval')
    pval_mat_sex$FDR    <- p.adjust(as.numeric(pval_mat_sex$pval), method = 'BH')
    print("Sex difference results:")
    print(pval_mat_sex)
    fwrite(pval_mat_sex, file = paste0(plot_dir, 'cs_organoid_nosub_sex_pval.csv'))
    
    # ASD difference
    pred_forstats <- pred[species=='human' & Sex=='Male']
    pred_means    <- pred_forstats[, .(mean_pred_time = median(pred_time)), by = .(DonorID, day, ASD_Ever_Dx, CellClass, batch, batchID)]
    pred_means    <- pred_means[,   .(mean_pred_time = median(mean_pred_time)), by = .(DonorID, day, ASD_Ever_Dx, CellClass, batchID)]
    pval_mat <- c()
    for (day_i in unique(pred_means$day)){
      for (ct_i in unique(pred_means$CellClass)){
        pval <- tryCatch(
          wilcox.test(mean_pred_time ~ ASD_Ever_Dx, data = pred_means[day==day_i & CellClass==ct_i])$p.value,
          error = function(e) NA_real_)
        pval_mat <- rbind(pval_mat, c(ct_i, day_i, pval))
      }
    }
    pval_mat_asd        <- as.data.table(pval_mat)
    names(pval_mat_asd) <- c('CellClass', 'day', 'pval')
    pval_mat_asd$FDR    <- p.adjust(as.numeric(pval_mat_asd$pval), method = 'BH')
    print("ASD difference results:")
    print(pval_mat_asd)
    fwrite(pval_mat_asd, file = paste0(plot_dir, 'cs_organoid_nosub_asd_pval.csv'))
  }
}
