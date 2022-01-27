Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# return the cosine distances for each changepoint in a trajectory
cosineDist <- function(traj) {
  bootstraps <- length(traj) / 4
  distances <- data.frame(cpPos1 = c(), dist = c())
  i <- 1
  for (k in 1:bootstraps) {
    if (!is.null(traj[[i+1]])) {
      for (j in traj[[i+1]]) {
        cosine_dist <- 1 - (as.numeric(lsa::cosine(traj[[i]][,j], traj[[i]][,j+1])))
        distances <- rbind(distances, data.frame(cpPos1=j, dist = cosine_dist))
      }
    }
    i <- i + 4
  }
  return (distances)
}

# find the amount by which each signature changes at all the changepoints in a trajectory
findSigChanges <- function(traj) {
  bootstraps <- length(traj) / 4
  sig_changes <- data.frame()
  i <- 1
  for (k in 1:bootstraps) {
    if (!is.null(traj[[i+1]])) {
      for (j in traj[[i+1]]) {
        changes <- as.data.frame(traj[[i]][,j+1]-traj[[i]][,j])
        changes_flip <- data.table::transpose(changes)
        colnames(changes_flip) <- rownames(changes)
        changes_flip$cpPos1 <- j
        sig_changes <- rbind(sig_changes, changes_flip)
      }
    }
    i <- i + 4
  }
  return (sig_changes)
}

# find the signatures that have the highest activity on either side of a changepoint in a trajectory
findDominantSigs <- function(traj) {
  bootstraps <- length(traj) / 4
  dominant_sigs <- data.frame()
  i <- 1
  for (k in 1:bootstraps) {
    if (!is.null(traj[[i+1]])) {
      for (j in traj[[i+1]]) {
        dominant_pre <- rownames(traj[[i]])[as.numeric(which(traj[[i]][,j] == max(traj[[i]][,j])))]
        dominant_post <- rownames(traj[[i]])[as.numeric(which(traj[[i]][,j+1] == max(traj[[i]][,j+1])))]
        sample <- data.frame(dominant_pre = dominant_pre, dominant_post = dominant_post, cpPos1 = j)
        dominant_sigs <- rbind(dominant_sigs, sample)
      }
    }
    i <- i + 4
  }
  return (dominant_sigs)
}

# make dataframe of information about all changepoints in a trajectory
summarizeChangepoints <- function(trajectory, sampleID) {
  cpPos <- assignChangepoints(trajectory,0)
  if (nrow(cpPos) > 0) {
    cpPos$start_chrom <- 1
    cpPos$end_chrom <- 1
    cpPos$start <- 1
    cpPos$end <- 1
    cpPos$mb_bin_start <- 1
    cpPos$mb_bin_end <- 1
    cpPos$sampleID <- sampleID
    for (i in 1:nrow(cpPos)) {
      if (!is.null(trajectory[['binData']]$actual_bin)) {
        cpPos$start[i] <- trajectory[['binData']]$start[trajectory[['binData']]$actual_bin==cpPos$cpPos1[i]]
        cpPos$end[i] <- trajectory[['binData']]$end[trajectory[['binData']]$actual_bin==cpPos$cpPos1[i]]
        cpPos$start_chrom[i] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$actual_bin==cpPos$cpPos1[i]]
        cpPos$end_chrom[i] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$actual_bin==cpPos$cpPos1[i]]
      }
      else {
        cpPos$start[i] <- trajectory[['binData']]$start[trajectory[['binData']]$bin==cpPos$cpPos1[i]]
        cpPos$end[i] <- trajectory[['binData']]$end[trajectory[['binData']]$bin==cpPos$cpPos1[i]]
        cpPos$start_chrom[i] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$bin==cpPos$cpPos1[i]]
        cpPos$end_chrom[i] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$bin==cpPos$cpPos1[i]]
      }
      mb_bin_start <- blank_genome %>%
        dplyr::filter(start_chrom == cpPos$start_chrom[i] & start == cpPos$start[i])
      cpPos$mb_bin_start[i] <- mb_bin_start$mb_bin
      mb_bin_end <- blank_genome %>%
        dplyr::filter(end_chrom == cpPos$end_chrom[i] & end == cpPos$end[i])
      cpPos$mb_bin_end[i] <- mb_bin_end$mb_bin
    }
    return (cpPos)
  }
  else {
    return (NULL)
  }
}

# make dataframe of information about all changepoints in a trajectory
summarizeChangepointsFreq <- function(trajectory, sampleID) {
  cpPos <- assignChangepoints(trajectory,0)
  if (nrow(cpPos) > 0) {
    cpPos$ccf_start <- 1
    cpPos$ccf_end <- 2
    cpPos$sampleID <- sampleID
    for (i in 1:nrow(cpPos)) {
      cpPos$ccf_start[i] <- colnames(trajectory[['mixtures']])[cpPos$cpPos1[i]]
      cpPos$ccf_end[i] <- colnames(trajectory[['mixtures']])[cpPos$cpPos2[i]]
    }
    return (cpPos)
  }
  else {
    return (NULL)
  }
}

buildChangepointDataFreq <- function(type) {
  cpPos <- data.frame()
  sig_changes_full <- data.frame()
  j <- 1

  for (i in c(1:330)) {
    if (file.exists(paste0("~/Desktop/CBSP2021/",type,"_vcf/",type,"_vcf",i,".Rdata"))) {
      print(i)

      load(paste0("~/Desktop/CBSP2021/",type,"_vcf/",type,"_vcf",i,".Rdata"))
      temp <- load(paste0("~/Desktop/CBSP2021/",type,"_vcf/",type,"_vcf",i,".Rdata"))
      melanoma_master <- c(melanoma_master, get(temp))

      # sig_changes <- findSigChanges(get(temp))
      # if (nrow(sig_changes)>0) {
      #   sig_changes$sampleID <- i
      #   sig_changes_full <- rbind(sig_changes_full, sig_changes)
      # }
      # cps <- summarizeChangepointsFreq(get(temp),i)
      # if (!is.null(cps)) {
      #   cps$tempID <- j
      #   cpPos <- rbind(cpPos, cps)
      # }
      # j <- j + 1
    }
  }
  if (nrow(cpPos) > 0) {
    changepoint_summaries <- cpPos %>%
      dplyr::rename(start_bin = cpPos1,
                    end_bin = cpPos2,
                    bootstrap_support = prob)
    changepoint_summaries <- cbind(changepoint_summaries, sig_changes_full)
    return (changepoint_summaries)
  }
  else {
    return (NULL)
  }
}

# build dataset of information for all changepoints for a cancer type
buildChangepointData <- function(type, blank_genome) {
  # initialize dataframes to combine
  cpPos <- data.frame()
  cosine_distances <- data.frame(cpPos1 = c(), dist = c(), sampleID = c(), tempID = c())
  sig_changes_full <- data.frame()
  dominant_sigs_full <- data.frame()
  mutation_counts <- c()
  #file_ids <- read.delim(paste0("~/Desktop/CBSP2021/",type,"_ids.txt"), header=FALSE)
  #colnames(file_ids) <- c("id")
  #cnas <- data.frame()

  j <- 1
  for (i in c(1:330)) {
    if (file.exists(paste0("~/Desktop/CBSP2021/",type,"/",type,i,".Rdata"))) {
      print(i)
      # temp_cna <- makeCnaDataset(trimws(as.character(file_ids$id[i]), which='both'), blank_genome)
      # temp_cna$sampleID <- i
      # cnas <- rbind(cnas, temp_cna)

      load(paste0("~/Desktop/CBSP2021/",type,"/",type,i,".Rdata"))
      temp <- load(paste0("~/Desktop/CBSP2021/",type,"/",type,i,".Rdata"))
      dominant_sigs <- findDominantSigs(get(temp))
      if (nrow(dominant_sigs)>0) {
        dominant_sigs$sampleID <- i
        dominant_sigs_full <- rbind(dominant_sigs_full, dominant_sigs)
      }
      sig_changes <- findSigChanges(get(temp))
      if (nrow(sig_changes)>0) {
        sig_changes$sampleID <- i
        sig_changes_full <- rbind(sig_changes_full, sig_changes)
      }
      dist <- cosineDist(get(temp))
      if (nrow(dist)>0) {
        dist$sampleID <- i
        dist$tempID <- j
        cosine_distances <- rbind(cosine_distances, dist)
      }
      cps <- summarizeChangepoints(get(temp),i)
      if (!is.null(cps)) {
        mutation_counts <- c(mutation_counts, rep(sum(colSums(get(temp)[['binData']][,7:102])), nrow(cps)))
        cps$tempID <- j
        cpPos <- rbind(cpPos, cps)
      }
      j <- j + 1
    }
  }
  if (nrow(cpPos)>0) {
    # find average cosine distances at each cp across bootstraps
    avg_distances <- cosine_distances %>%
      dplyr::group_by(sampleID, cpPos1) %>%
      dplyr::summarize(avg_dist = mean(dist))
    # find average signature exposure changes at each cp across bootstraps
    avg_sig_changes <- sig_changes_full %>%
      dplyr::group_by(sampleID, cpPos1) %>%
      dplyr::summarise_at(dplyr::vars(colnames(sig_changes_full)[1]:colnames(sig_changes_full)[ncol(sig_changes_full)-2]), mean) %>%
      dplyr::select(-c(sampleID,cpPos1))
    # find average of dominant signatures at each cp across bootstraps
    avg_dominant_sigs <- dominant_sigs_full %>%
      dplyr::group_by(sampleID, cpPos1) %>%
      dplyr::summarize_at(dplyr::vars(dominant_pre:dominant_post), Mode)

    changepoint_summaries <- cpPos %>%
      dplyr::rename(start_bin = cpPos1,
                    end_bin = cpPos2,
                    bootstrap_support = prob) %>%
      dplyr::mutate("avg_cosine_dist" = avg_distances$avg_dist,
                    "dominant_pre" = avg_dominant_sigs$dominant_pre,
                    "dominant_post" = avg_dominant_sigs$dominant_post)
    changepoint_summaries <- cbind(changepoint_summaries, avg_sig_changes[,2:length(colnames(avg_sig_changes))])
    return (changepoint_summaries)
  }
  else {
    return (NULL)
  }


}

# KDE
changepoint_positions <- c()
for (i in 1:nrow(melanoma_changepoints)) {
  changepoint_positions <- c(changepoint_positions, seq(from=melanoma_changepoints$mb_bin_start[i], to=melanoma_changepoints$mb_bin_end[i]))
}
changepoint_positions <- as.data.frame(changepoint_positions)
kde <- ks::kde(changepoint_positions, h=4)
kde_data <- data.frame("eval" = kde$eval.points, "estimate" = kde$estimate)
ggplot2::ggplot() +
  ggplot2::geom_histogram(data = changepoint_positions, mapping=ggplot2::aes(x=changepoint_positions), binwidth=2, fill='orange') +
  ggplot2::geom_line(data = kde_data, mapping=ggplot2::aes(x = eval, y = estimate*2500)) +
  ggplot2::theme_bw() +
  ggplot2::scale_x_continuous(breaks = change_bins, labels = chr_labels) +
  ggplot2::labs(x = "Chromosome", y = "Changepoint Frequency", fill = 'Distribution') +
  ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank())

buildOverlapsData <- function(cpPos, bandwidth) {
  changepoints <- c()
  overlaps_full <- data.frame('sampleID' = c(),'sd_min' = c(), 'sd_max' = c())
  for (i in 1:nrow(cpPos)) {
    cp_range <- rep(cpPos$mb_bin_start[i]:cpPos$mb_bin_end[i])
    if (length(cp_range)>1) {
      kde <- ks::kde(cp_range, bandwidth)
      overlaps <- data.frame(sampleID = cpPos$sampleID[i],
                             sd_min = kde$eval.points[which.min(abs(kde$estimate-max(kde$estimate)))] - sd(kde$eval.points),
                             sd_max =  kde$eval.points[which.max(abs(kde$estimate+max(kde$estimate)))] + +sd(kde$eval.points))
      overlaps_full <- rbind(overlaps_full, overlaps)
    }
    else {
      cp_range <- c(cpPos$mb_bin_start[i], (cpPos$mb_bin_start[i]+.5), (cpPos$mb_bin_start[i]+.99))
      kde <- ks::kde(cp_range, bandwidth)
      overlaps <- data.frame(sampleID = cpPos$sampleID[i],
                             sd_min = kde$eval.points[which.min(abs(kde$estimate-max(kde$estimate)))] - sd(kde$eval.points),
                             sd_max =  kde$eval.points[which.max(abs(kde$estimate+max(kde$estimate)))] + +sd(kde$eval.points))
      overlaps_full <- rbind(overlaps_full, overlaps)
    }

  }
  overlaps_full$tempID <- cpPos$tempID

  return (overlaps_full)
}




# PLOTTING
plotChangepointSummary <- function(cpPos, overlaps_full=NULL, blank_genome, title, subtitle) {

  change_bins <- c(1)
  blank_genome <- readr::read_csv("~/Desktop/CBSP2021/Thy-AdenoCA_pooled.csv")[,1:4]
  blank_genome$mb_bin <- rep(1:nrow(blank_genome))
  for (i in 2:nrow(blank_genome)-1){
    if (blank_genome$start_chrom[i] < blank_genome$start_chrom[i+1]) {
      change_bins <- c(change_bins, blank_genome$mb_bin[i+1])
    }
  }
  chr_labels <- as.character(c(1:22, "X", "Y"))

  g <- ggplot2::ggplot(data = blank_genome, ggplot2::aes(x = mb_bin)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = change_bins, labels = chr_labels) +
    ggplot2::labs(x = "Chromosome", title = title,
                  subtitle = subtitle, y = "Sample") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank())

  for (i in 1:length(change_bins)-1) {
    if (i %% 2 != 0) {
      g <- g + ggplot2::annotate("rect", xmin=change_bins[i], xmax = change_bins[i+1],
                                 ymin=-Inf, ymax=Inf, alpha=0.3, fill='lightgrey')
    }
  }
  if (!is.null(overlaps_full)) {
    for (i in 1:nrow(overlaps_full)) {
      g <- g + ggplot2::annotate("rect", xmin = overlaps_full$sd_min[i], xmax = overlaps_full$sd_max[i],
                                 ymin=overlaps_full$tempID[i], ymax=overlaps_full$tempID[i]+1, fill="orange", alpha=0.3)
    }
  }
  for (i in 1:nrow(cpPos)) {
    g <- g + ggplot2::annotate("rect", xmax=cpPos$mb_bin_end[i], xmin=cpPos$mb_bin_start[i],
                               ymin=cpPos$tempID[i], ymax=cpPos$tempID[i]+1, alpha=cpPos$bootstrap_support[i], fill = "red")
  }
  g <- g + ggplot2::annotate("rect", xmax=3100, xmin=5,
                             ymin=max(cpPos$tempID), ymax = max(cpPos$tempID)+1, alpha=0, fill='white')
  g
}

# IDENTIFY CHANGEPOINT PEAKS

calculatePeaks <- function(overlaps) {
  peaks <- data.frame("start" = c(), "stop" = c(), "n_cps" = c())
  for (k in c(1:3109)) {
    peak <- nrow(overlaps %>%
                   dplyr::filter((sd_min>=k & sd_min<=(k+5)) | (sd_max>=k & sd_max<=(k+5)) | (sd_min<=k & sd_max >= (k+5))))
    peaks <- rbind(peaks, data.frame("start" = k, "stop" = k+5, "n_cps" = peak))
  }
  peaks <- peaks %>%
    dplyr::filter(n_cps > 2)

  return (peaks)
}

wranglePeaks <- function(peaks, overlaps, changepoints, n, y) {
  peaks <- peaks %>%
    dplyr::mutate('end' = start+5,
                  'sample_count' = 0)
  overlaps$peak <- 0
  for (k in 1:nrow(peaks)) {
    peak <- overlaps %>%
      dplyr::filter((sd_min>=peaks$start[k] & sd_min<=peaks$end[k]) | (sd_max>=peaks$start[k] & sd_max<=peaks$end[k]) | (sd_min<=peaks$start[k] & sd_max >=peaks$end[k]))
    sample_count <- peak %>%
      dplyr::group_by(sampleID) %>%
      dplyr::summarize(n = dplyr::n_distinct())
    peaks$sample_count[k] <- nrow(sample_count)
    for (j in 1:nrow(overlaps)) {
      if (overlaps$sd_min[j] %in% peak$sd_min) {
        overlaps$peak[j] <- peaks$start[k]
      }
    }
  }
  sig_changes <- changepoints[,15:ncol(changepoints)]
  max_increase <- names(sig_changes)[apply(sig_changes, 1, which.max)]
  max_decrease <- names(sig_changes)[apply(sig_changes, 1, which.min)]
  sig_changes <- sig_changes %>%
    dplyr::mutate(max_increase = max_increase,
                  max_decrease = max_decrease)

  plot_table <- changepoints %>%
    cbind(overlaps[,c(2,3,5)]) %>%
    cbind(sig_changes[,c(ncol(sig_changes),ncol(sig_changes)-1)])
  plot_table <- plot_table %>%
    dplyr::group_by(peak) %>%
    dplyr::summarize(sample_support = dplyr::n_distinct(sampleID),
                     max_increase = Mode(max_increase),
                     max_decrease = Mode(max_decrease),
                     mean_dist = mean(avg_cosine_dist)) %>%
    dplyr::filter(peak != 0,
                  sample_support > 2) %>%
    dplyr::mutate(y = y)
  for (sig in colnames(TrackSig::cosmicV3)) {
    plot_table <- plot_table %>%
      dplyr::mutate(!!sig := dplyr::case_when((max_increase == sig & max_decrease == sig) ~ 1,
                                              (max_decrease == sig & max_increase != sig) ~ .5,
                                              (max_increase == sig & max_decrease != sig) ~ .5,
                                              (max_increase != sig & max_decrease != sig) ~ 0))
  }

  # plot_table <- plot_table %>%
  #   tidyr::pivot_longer(cols=colnames(plot_table)[7:84], names_to="sig", values_to="value") %>%
  #   dplyr::filter(value > 0) %>%
  #   dplyr::arrange(dplyr::desc(value))
  #
  # change_order <- plot_table %>%
  #   dplyr::filter(value < 1)
  #
  # plot_table <- plot_table %>%
  #   dplyr::filter(value > 0.5)
  #
  # for (i in seq(1, (nrow(change_order)-1), 2)) {
  #   change_order$sig[i] <- change_order$max_increase[i]
  #   change_order$sig[i+1] <- change_order$max_decrease[i]
  # }
  #
  # plot_table <- rbind(plot_table, change_order)

  return (plot_table)
}

buildPlottingData <- function(type, blank_genome, peaks, n, y) {
  changepoints <- buildChangepointData(type, blank_genome)
  overlaps <- buildOverlapsData(changepoints, 4)
  plot_table <- wranglePeaks(peaks, overlaps, changepoints, n, y)
  return (plot_table)
}



type_breaks <- c(1,150,300,450,600,750,900,1050,1200,
                 1350,1500,1650)
cancer_types <- c("Melanoma (70/107)", "Lung-SCC (25/48)", "Eso-AdenoCA (45/97)",
                  "Colorect-AdenoCA (21/60)", "Bladder-TCC (13/23)", "Stomach-AdenoCA (6/67)",
                    "Lymph-BNHL (25/106)",
                  "Uterus-AdenoCA (8/51)", "CNS-GBM (21/41)",
                  "Bone-Osteosarc (6/39)",
                  "Lymph-CLL (58/95)", "Cervix (3/20)")

changepoints <- buildChangepointData('cervix', blank_genome)
peaks <- cervix_peaks
overlaps <- buildOverlapsData(changepoints, 4)
peaks <- peaks %>%
  dplyr::mutate('end' = start+5,
                'sample_count' = 0)
overlaps$peak <- 0
for (k in 1:nrow(peaks)) {
  peak <- overlaps %>%
    dplyr::filter((sd_min>=peaks$start[k] & sd_min<=peaks$end[k]) | (sd_max>=peaks$start[k] & sd_max<=peaks$end[k]) | (sd_min<=peaks$start[k] & sd_max >=peaks$end[k]))
  sample_count <- peak %>%
    dplyr::group_by(sampleID) %>%
    dplyr::summarize(n = dplyr::n_distinct())
  peaks$sample_count[k] <- nrow(sample_count)
  for (j in 1:nrow(overlaps)) {
    if (overlaps$sd_min[j] %in% peak$sd_min) {
      overlaps$peak[j] <- peaks$start[k]
    }
  }
}
plot <- wranglePeaks(peaks, overlaps, changepoints, 20, 1)
changepoints$peak <- overlaps$peak
changepoints %>%
  dplyr::filter(peak %in% plot$peak) %>%
  dplyr::summarize(n = dplyr::n_distinct(sampleID))


combined_plot <- combined_plot[!is.na(combined_plot$y), ]
combined_plot <- combined_plot %>%
  dplyr::mutate(y = dplyr::case_when(y == 600 ~ 600,
                                     y == 1 ~ 1,
                                     y == 150 ~ 150,
                                     y == 300 ~ 300,
                                     y == 450 ~ 450,
                                     y == 750 ~ 750,
                                     y == 900 ~ 900,
                                     y == 1050 ~ 1050,
                                     y == 1200 ~ 1200,
                                     y == 1650 ~ 1350,
                                     y == 1950 ~ 1500,
                                     y == 2100 ~ 1650))
combined_plot <- rbind(melanoma_plot, scc_plot, eso_plot, colorect_plot, bladder_plot, stomach_plot,
                       bnhl_plot, uterus_plot, gbm_plot, breast_plot,
                       panc_plot, bone_plot, prostate_plot, cll_plot, cervix_plot)
combined_plot <- combined_plot %>%
  dplyr::filter(sample_support >= 3)
combined_plot$sig <- factor(combined_plot$sig, levels=colnames(TrackSig::cosmicV3))

default_palette <- ggpubr::get_palette(palette="igv", k=24)
default_palette <- default_palette[-c(21,22)]

changepoint_summary <- ggplot2::ggplot(data = combined_plot, ggplot2::aes(x = peak, y = y)) +
  scatterpie::geom_scatterpie(ggplot2::aes(x=peak, y=y, r=combined_plot$sample_support*1.5),
                              pie_scale=1, color=NA, data=combined_plot, cols="sig", long_format = TRUE
                              ) +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  ggplot2::scale_fill_manual(values=default_palette) +
  #ggsci::scale_fill_igv("default") +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol=1)) +
  ggplot2::theme(legend.position='right',
                 panel.grid.minor.y = ggplot2::element_blank(),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(size=8),
                 axis.text.y = ggplot2::element_text(size=10),
                 axis.title.x = ggplot2::element_text(size=12),
                 axis.title.y = ggplot2::element_text(size=12),
                 legend.title = ggplot2::element_text(size=12),
                 legend.text = ggplot2::element_text(size=10),
                 plot.margin = ggplot2::margin(0,0,0,0)) +
  ggplot2::scale_x_continuous(breaks = change_bins, labels = chr_labels) +
  ggplot2::scale_y_continuous(breaks=type_breaks, labels=cancer_types) +
  ggplot2::labs(x = "Chromosome", fill="Signature", y = "Cancer Type")

ggplot2::ggsave(filename="updated_changepoint_summary.pdf", plot=ggplot2::last_plot(), device="pdf", path="~/Desktop/CBSP2021/plots",
                scale=1, width = 35, height = 18, units = "cm", dpi=300)

melanoma_peaks <- data.frame('start' = c(47,191,481,565,703,898,908,1008,1050,1113,1195,1243,
                                         1389,1521,1531,1167,1760,1797,1881,1900,1928,2055,2141,2198,2294,
                                         2496,2560,2696,2826))
lung_scc_peaks <- data.frame('start' = c(47,539,671,689,913,1109,1561,2114,2149,2180,2188,2214,2241,2873,2885,2897))
bladder_peaks <- data.frame('start' = c(47,539,671,689,913,1109,1561,2149,2162,2180,2188,2214,2241,2672,2788,2873,2885,2897,2928))
eso_peaks <- data.frame('start' = c(42,256,674,689,707,1532,1562,1576,1620,2051,2084,2092,2138,2181,2189,2256,2290,2348,2352,2358,
                                    2488,2495,2563,2581,2662,2790,2824,2836,2876,2895,2922,2938,2946))
lymph_bnhl_peaks <- data.frame('start' = c(707,911,2181,2194,2201,2208,2248,2275,2290,2308,2840,2864,2906))
prostate_peaks <- data.frame('start' = c(1392,1703))
stomach_peaks <- data.frame('start' = c(47,2292,2907,2917))
cll_peaks <- data.frame('start' = c(426,457,485,524,580,601,612,624,656,668,675,832,846,861,867,875,908,1097,1367,
                                    1435,1514,1536,1557,1581,1607,1776,2182,2238,2264,2331,2494,2623,2657,2734,2835))
cns_peaks <- data.frame('start' = c(66,251,704,905,1170,1290,1351,1875,1926,2142,2185,2249,2269,2284,2335,
                                    2855,2913))
colorect_peaks <- data.frame('start'= c(42,58,104,191,388,478,567,596,604,702,754,826,897,908,
                                        1054,1110,1193,1214,1237,1249,1293,1391,1452,1532,1669,1725,
                                        1807,1899,1931,1967,2062,2142,2179,2198,2254,2294,2480,2615,
                                        2659,2696,2830,2967))
uterus_peaks <- data.frame('start' = c(47,191,481,705,891,1052,1110,1195,1243,1405,1536,1669,1736,1818,1895,2015,
                                       2062,2145,2188,2297,2419,2480,2695,2782,2824))
cervix_peaks <- data.frame('start' = c(2933))
breast_peaks <- data.frame('start' = c(1205,2894))
panc_adenoca_peaks <- data.frame('start' = c(2022,2494,2536,2664))
bone_osteosarc_peaks <- data.frame('start' =c(1088,1102,2068))

peaks <- melanoma_peaks
changepoints <- melanoma_changepoints
table <- overlaps

for (i in 1:nrow(peaks)) {
  data <- changepoints %>%
    dplyr::mutate("peak" = table$peak) %>%
    dplyr::filter(peak == 47) %>%
    dplyr::select("SBS1":"peak", "sampleID", "mb_bin_start")
  data <- data %>%
    dplyr::mutate('id' = rep(1:nrow(data)))
  data <- data %>%
    base::subset(select=c((ncol(data)-3),(ncol(data)-2),(ncol(data)-1),ncol(data),1:(ncol(data)-4)))

  data <- data %>%
    tidyr::pivot_longer(cols=colnames(data)[5:length(colnames(data))],
                        values_to = "activity_change", names_to="sig")
  data <- data %>%
    dplyr::mutate(direction = dplyr::if_else(activity_change<0, "decrease", "increase"),
                  alpha = dplyr::if_else((abs(activity_change) >=0.01), 1, 0),
                  sig = factor(sig, levels=c(colnames(changepoints)[15:length(colnames(changepoints))])))
  data$direction <- factor(data$direction, levels = c("decrease", "increase"))

  ggplot2::ggplot(data = data, ggplot2::aes(x=direction, y=abs(activity_change)*100, fill=sig)) +
    ggplot2::geom_bar(position='dodge', width=1, stat='identity') +
    ggplot2::scale_fill_manual(values=melanoma_palette) +
    #ggsci::scale_fill_igv(palette='default') +
    ggplot2::scale_x_discrete(limits=rev(levels(data$direction))) +
    ggplot2::coord_polar() +
    ggplot2::labs(y = "Exposure change (%)", fill = "Signature", x = data$sampleID) +
    ggplot2::facet_wrap(.~sampleID, ncol=10) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size=11, angle=90),
                   axis.text.y = ggplot2::element_text(size=10),
                   legend.title = ggplot2::element_text(size=16),
                   legend.text = ggplot2::element_text(size=14),
                   axis.title.y = ggplot2::element_text(size=16),
                   axis.title.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=22),
                   strip.background = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_blank())
                   #strip.text.x = ggplot2::element_text(size=16))

  ggplot2::ggsave(filename=paste0("melanoma47_paper.pdf"), path="~/Desktop/CBSP2021/plots", device="pdf",
         plot=ggplot2::last_plot(), scale=1, width=45, height=25, unit="cm")
}

missing_data <- data.frame("type" = types, "med_cps" = med_cps, "med_dist" = med_dists)
med_cps <- c()
med_dists <- c()
types <- c('bladder', 'bone-osteosarc', 'cervix',
           'lymph-cll', 'cns-gbm', 'colorect', 'eso-adenoca',
           'lung-scc',
           'lymph-bnhl', 'melanoma', 'stomach',
           'uterus-adenoca')
peaks <- list(bladder_peaks, bone_osteosarc_peaks, cervix_peaks,
           cll_peaks, cns_peaks, colorect_peaks, eso_peaks, lung_scc_peaks,
           lymph_bnhl_peaks, melanoma_peaks, stomach_peaks,
           uterus_peaks)
ns <- c(23, 39, 20, 95, 41, 60, 97, 48, 106, 107, 67, 51)

for (j in 1:length(types)) {
  # cps <- buildChangepointData(t, blank_genome)
  # med_ncps <- cps %>%
  #   dplyr::group_by(sampleID) %>%
  #   dplyr::summarize(n = dplyr::n())
  # med_cps <- c(med_cps, median(med_ncps$n))
  # med_dists <- c(med_dists, median(cps$avg_cosine_dist))
  peaks <- as.data.frame(peaks[[j]])
  colnames(peaks) <- c('start')

  uterus_plot <- buildPlottingData(types[j], blank_genome, uterus_peaks, ns[j], type_breaks[j])
  revised_plot <- rbind(revised_plot, plot)
}

combined_plot <- rbind(melanoma_plot, scc_plot, eso_plot, colorect_plot,
                      bladder_plot, stomach_plot, bnhl_plot, uterus_plot,
                      cns_plot, bone_plot, cll_plot, cervix_plot)
readr::write_csv(combined_plot, "revised_plot.csv")
