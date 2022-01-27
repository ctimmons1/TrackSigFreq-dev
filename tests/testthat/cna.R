mbBinMaster <- function(master, binSize) {
  start_chrom <- NULL
  master$mb_bin <- blank_genome$mb_bin
  # initialize empty dataframe
  counts <- data.frame()
  for (i in c(1:23)) {
    # find mutation counts in each row
    # bin sex chromosomes together
    if (i == 23) {
      master_subset <- master %>%
        dplyr::filter(start_chrom>=i)
      master_subset <- master_subset %>%
        dplyr::mutate(rowsum = rowSums(master_subset[,6:101]),
                      bin = 1) %>%
        base::subset(select = c(1:5,102,103,104,6:101))

      bin <- 1
      sums <- 0
      # add up counts in each row until desired bin size is reached
      for (i in 1:nrow(master_subset)) { # nolint
        sums <- sums + master_subset$rowsum[i]
        master_subset$bin[i] <- bin
        if (sums >= binSize) {
          sums <- 0
          bin <- bin + 1
        }
      }
      counts <- rbind(counts, master_subset)
    }
    else {
      # bin autosomal chromosomes individually
      master_subset <- master %>%
        dplyr::filter(start_chrom==i)
      master_subset <- master_subset %>%
        dplyr::mutate(rowsum = rowSums(master_subset[,6:101]),
                      bin = 1) %>%
        base::subset(select = c(1:5,102,103,104,6:101))

      bin <- 1
      sums <- 0
      # add up counts in each row until desired bin size is reached
      for (i in 1:nrow(master_subset)) {
        sums <- sums + master_subset$rowsum[i]
        master_subset$bin[i] <- bin
        if (sums >= binSize) {
          sums <- 0
          bin <- bin + 1
        }
      }
      counts <- rbind(counts, master_subset)
    }
  }
  # define bins in relation to the whole genome
  for (i in c(2:23)) {
    if (i == 23) {
      counts$bin[counts$start_chrom>=i] <- counts$bin[counts$start_chrom>=i] +
        max(counts$bin[counts$start_chrom==i-1])
    }
    else {
      counts$bin[counts$start_chrom==i] <- counts$bin[counts$start_chrom==i] +
        max(counts$bin[counts$start_chrom==i-1])
    }
  }

  return (counts)
}


makeCnaDataset <- function(fileID, blank_genome) {
  cna_dataset <- data.frame()

    file <- read.delim(paste0("~/Desktop/CBSP2021/consensus.20170119.somatic.cna.annotated/",fileID,".consensus.20170119.somatic.cna.annotated.txt"))
    file <- file %>%
      dplyr::select(chromosome, start, end, total_cn) %>%
      dplyr::mutate(mb_bin_start = 1,
                    mb_bin_end = 1,
                    chromosome = dplyr::case_when(chromosome == "X" ~ 23,
                                                  chromosome == "Y" ~ 24,
                                                  TRUE ~ as.numeric(chromosome)),
                    total_cn = dplyr::case_when(total_cn == NA ~ 2,
                                                 TRUE ~ as.numeric(total_cn)))

    for (i in c(1:24)) {
      blank_subset <- blank_genome %>%
        dplyr::filter(start_chrom==i)
      file_subset <- file %>%
        dplyr::filter(chromosome==i)
      if (nrow(file_subset)>0) {
        for (j in 1:nrow(file_subset)) {
          file_subset$mb_bin_start[j] <- blank_subset$mb_bin[(floor(file_subset$start[j])/1e6)+1] + ((file_subset$start[j] - (trunc((file_subset$start[j]/1e6),digits=0)*1e6))/1e6)
          file_subset$mb_bin_end[j] <- blank_subset$mb_bin[(ceiling(file_subset$end[j])/1e6)+1] + ((file_subset$end[j] - (trunc((file_subset$end[j]/1e6),digits=0)*1e6))/1e6)
        }
        cna_dataset <- rbind(cna_dataset, file_subset)
      }
    }
    cna_dataset <- cna_dataset %>%
      tidyr::drop_na()
    return (cna_dataset)

}

# Plot changepoints and CNAs
colnames(file_ids) <- c("id")
cnas <- makeCnaDataset(file_ids$id[1], blank_genome)
cpPos <- summarizeChangepoints(traj,18)

g <- ggplot2::ggplot(data = cnas, ggplot2::aes(x = mb_bin)) +
  ggplot2::theme_bw() +
  ggplot2::scale_x_continuous(breaks = change_bins, labels = chr_labels) +
  ggplot2::labs(x = "Chromosome", y = "Copy Number", title="Eso-AdenoCA 19 (pval=0.051)") +
  ggplot2::theme(
                 panel.grid.major.y = ggplot2::element_blank(),
                 panel.grid.minor.y = ggplot2::element_blank(),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 panel.grid.major.x = ggplot2::element_blank(),
                 legend.position = 'none')
for (i in 1:length(change_bins)-1) {
  if (i %% 2 != 0) {
    g <- g + ggplot2::annotate("rect", xmin=change_bins[i], xmax = change_bins[i+1],
                               ymin=-Inf, ymax=Inf, alpha=0.3, fill='lightgrey')
  }
}
for (i in 1:nrow(cpPos)) {
  g <- g + ggplot2::annotate("rect", xmax=cpPos$mb_bin_end[i], xmin=cpPos$mb_bin_start[i], ymin=-Inf, ymax=Inf, alpha=cpPos$prob[i], fill = 'red')
}
g <- g + ggplot2::annotate("rect", xmax=250, xmin=1,
                           ymin=max(cpPos$sampleID), ymax = max(cpPos$sampleID)+1, alpha=0, fill='white')
g <- g + ggplot2::geom_step(ggplot2::aes_string(x = cnas$mb_bin_start, y = cnas$total_cn, alpha=0.7)) +
  ggplot2::scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10), labels = as.character(rep(1:10)), limits = c(0,10))

g

