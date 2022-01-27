generateSampleOverlaps <- function(fileID, traj, blank_genome) {
  # get changepoints
  changepoints <- summarizeChangepoints(traj, 1) %>%
    dplyr::mutate(width = mb_bin_end - mb_bin_start)
  # get CNAs
  cnas <- makeCnaDataset(fileID, blank_genome) %>%
    dplyr::mutate_at(c("total_cn"), ~replace(., is.na(.), 2)) %>%
    dplyr::filter(total_cn != 2) %>%
    tidyr::drop_na()

  if (nrow(cnas) > 0) {
    # vector to store overlaps of all changepoints in the sample
    sample_overlaps <- c()
    # find the proportion of overlap between each changepoint and its closest CNA
    for (i in 1:nrow(changepoints)) {
      # find closest CNA on either side of changepoint
      sample_cna_lower <- which.min(abs(changepoints$mb_bin_start[i]-cnas$mb_bin_end))
      sample_cna_higher <- which.min(abs((changepoints$mb_bin_end[i])-cnas$mb_bin_start))
      # find which CNA on either side of the changepoint is closest to the changepoint
      if (abs((cnas$mb_bin_start[sample_cna_higher]-changepoints$mb_bin_end[i])) < abs((cnas$mb_bin_end[sample_cna_lower]-changepoints$mb_bin_start[i]))) {
        sample_cna_index <- sample_cna_higher
      } else {
        sample_cna_index <- sample_cna_lower
      }
      # range of coordinates containing the changepoint and the CNA
      sample_cna_range <- IRanges::IRanges(rep(cnas$mb_bin_start[sample_cna_index]:cnas$mb_bin_end[sample_cna_index]))
      sample_cp_range <- IRanges::IRanges(rep(changepoints$mb_bin_start[i]:changepoints$mb_bin_end[i]))
      cp_overlaps <- IRanges::overlapsAny(sample_cna_range, sample_cp_range)
      if (TRUE %in% cp_overlaps) {
        cp_overlaps <- 1
      }
      # proportions of overlaps for each changepoint in the sample
      sample_overlaps <- c(sample_overlaps, cp_overlaps)
    }
    # we have the overlaps for all the changepoints in the sample
    # find the proportion of changepoints in the sample that overlap with a CNA
    sample_prop <- length(sample_overlaps[sample_overlaps>0])/length(sample_overlaps)
    return (sample_prop)
  }
  else {
    return (rep(0, each=nrow(changepoints)))
  }

}

generateRandomOverlaps <- function(fileID, traj, blank_genome) {
  # get changepoints -- to know number of cps in the sample
  changepoints <- summarizeChangepoints(traj, 1) %>%
    dplyr::mutate(width = mb_bin_end - mb_bin_start)
  # get cnas
  cnas <- makeCnaDataset(fileID, blank_genome) %>%
    dplyr::mutate_at(c("total_cn"), ~replace(., is.na(.), 2)) %>%
    dplyr::filter(total_cn != 2) %>%
    tidyr::drop_na()

  if (nrow(cnas) > 0) {
    null <- c()
    for (j in c(1:10000)) {
      # randomize positions of changepoints
      changepoints$mb_bin_start <- sample(rep(1:3113), size=nrow(changepoints), replace=FALSE)
      changepoints$mb_bin_end <- changepoints$mb_bin_start + changepoints$width
      # initialize vector to store overlaps for all changepoints in random sample
      random_overlaps <- c()
      for (i in 1:nrow(changepoints)) {
        # find closest CNA on either side of random changepoint
        sample_cna_lower <- which.min(abs(changepoints$mb_bin_start[i]-cnas$mb_bin_end))
        sample_cna_higher <- which.min(abs((changepoints$mb_bin_end[i])-cnas$mb_bin_start))
        # find which CNA is closest to the changepoint and use that to calculate overlap
        if (abs((cnas$mb_bin_start[sample_cna_higher]-changepoints$mb_bin_end[i])) < abs((cnas$mb_bin_end[sample_cna_lower]-changepoints$mb_bin_start[i]))) {
          sample_cna_index <- sample_cna_higher
        } else {
          sample_cna_index <- sample_cna_lower
        }
        # ranges of coordinates containing the random changepoint and its closest CNA
        sample_cna_range <- rep(cnas$mb_bin_start[sample_cna_index]:cnas$mb_bin_end[sample_cna_index])
        sample_cp_range <- rep(changepoints$mb_bin_start[i]:changepoints$mb_bin_end[i])

        sample_cna_range <- IRanges::IRanges(rep(cnas$mb_bin_start[sample_cna_index]:cnas$mb_bin_end[sample_cna_index]))
        sample_cp_range <- IRanges::IRanges(rep(changepoints$mb_bin_start[i]:changepoints$mb_bin_end[i]))
        cp_overlaps <- IRanges::overlapsAny(sample_cna_range, sample_cp_range)
        if (TRUE %in% cp_overlaps) {
          cp_overlaps <- 1
        }
        random_overlaps <- c(random_overlaps, cp_overlaps)
      }
      # proportion of changepoints that overlap in the random sample
      random_prop <- length(random_overlaps[random_overlaps>0])/length(random_overlaps)
      # add this to null distribution of 10k samples
      null <- c(null, random_prop)
    }
    return (null)
  }

  else {
    return (rep(0, each=10000))
  }


}

testOverlapSignificance <- function(random_overlaps, sample_overlap) {
  p <- length(random_overlaps[random_overlaps>sample_overlap])/length(random_overlaps)
  return (p)
}

overlapTest <- function(fileID, traj, blank_genome) {
  random_overlaps <- generateRandomOverlaps(fileID, traj, blank_genome)
  sample_overlaps <- generateSampleOverlaps(fileID, traj, blank_genome)
  pval <- testOverlapSignificance(random_overlaps, sample_overlaps)
  return (pval)
}

testAllSamples <- function(type, blank_genome) {
  pvals_all <- c()
  ids <- c()


  for (i in c(1:330)) {
    if (file.exists(paste0("~/Desktop/CBSP2021/",type,"/",type,i,".Rdata"))) {
      if (file.exists(paste0("~/Desktop/CBSP2021/",type,"_ids.txt"))) {
        print(i)
        file_ids <- read.delim(paste0("~/Desktop/CBSP2021/",type,"_ids.txt"), header=FALSE)
        colnames(file_ids) <- c("id")

        load(paste0("~/Desktop/CBSP2021/",type,"/",type,i,".Rdata"))
        temp <- load(paste0("~/Desktop/CBSP2021/",type,"/",type,i,".Rdata"))

        # check for changepoints
        if (!is.null(get(temp)[[2]])||!is.null(get(temp)[[6]])||!is.null(get(temp)[[10]])||!is.null(get(temp)[[14]])||!is.null(get(temp)[[18]])) {
          # ||!is.null(get(temp)[[22]])||!is.null(get(temp)[[26]])||!is.null(get(temp)[[30]])||!is.null(get(temp)[[34]])||!is.null(get(temp)[[38]])) {
          fileID <- trimws(as.character(file_ids$id[i]), which='both')
          if (file.exists(paste0("~/Desktop/CBSP2021/consensus.20170119.somatic.cna.annotated/",fileID,".consensus.20170119.somatic.cna.annotated.txt"))) {
            # test for overlap
            ids <- c(ids, i)
            pval <- overlapTest(fileID, get(temp), blank_genome)
            pvals_all <- c(pvals_all, pval)
          }
        }
      }
    }
  }

  results <- data.frame('id' = ids, 'pval' = pvals_all)
  return (results)
}
