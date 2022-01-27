findKatDistances <- function(type, abbr) {
  kataegis_full <- readr::read_delim("~/Desktop/CBSP2021/combined_results_03032017.full.summary",
                                     delim="\t", col_names=TRUE)
  kataegis_full <- kataegis_full %>%
    tidyr::separate(sample, into=c("type", "id"), sep="_", remove=TRUE)

  kat_type <- kataegis_full %>%
    dplyr::filter(type %in% abbr) %>%
    dplyr::mutate(mb_bin_start = 1,
                  mb_bin_end = 1)

  kat_subset <- data.frame()
  for (i in c(1:24)) {
    blank_subset <- blank_genome %>%
      dplyr::filter(start_chrom==i)
    file_subset <- kat_type %>%
      dplyr::filter(chr==i)
    if (nrow(file_subset)>0) {
      for (j in 1:nrow(file_subset)) {
        file_subset$mb_bin_start[j] <- blank_subset$mb_bin[(floor(file_subset$start[j])/1e6)+1] + ((file_subset$start[j] - (trunc((file_subset$start[j]/1e6),digits=0)*1e6))/1e6)
        file_subset$mb_bin_end[j] <- blank_subset$mb_bin[(ceiling(file_subset$end[j])/1e6)+1] + ((file_subset$end[j] - (trunc((file_subset$end[j]/1e6),digits=0)*1e6))/1e6)
      }
      kat_subset <- rbind(kat_subset, file_subset)
    }
  }

  # make changepoint dataset
  changepoints <- buildChangepointData(type, blank_genome)
  # read in ids
  file_ids <- read.delim(paste0("~/Desktop/CBSP2021/",type,"_ids.txt"), header=FALSE)
  colnames(file_ids) <- c("id")
  file_ids$id <- trimws(as.character(file_ids$id), which='both')


  changepoints$id <- file_ids$id[changepoints$sampleID]
  distances <- c()
  for (i in 1:nrow(changepoints)) {
    sample_kat <- kat_subset %>%
      dplyr::filter(id == changepoints$id[i])
    if (nrow(sample_kat) > 0) {
      # find closest kataegis event to changepoint
      kat_before <- which.min(abs(changepoints$mb_bin_start[i]-sample_kat$mb_bin_end))
      kat_after <- which.min(abs((changepoints$mb_bin_end[i])-sample_kat$mb_bin_start))
      # find which CNA on either side of the changepoint is closest to the changepoint
      if (abs((sample_kat$mb_bin_start[kat_after]-changepoints$mb_bin_end[i])) < abs((sample_kat$mb_bin_end[kat_before]-changepoints$mb_bin_start[i]))) {
        closest_kat <- kat_after
      } else {
        closest_kat <- kat_before
      }
      kat_range <- IRanges::IRanges(rep(sample_kat$mb_bin_start[closest_kat]:sample_kat$mb_bin_end[closest_kat]))
      cp_range <- IRanges::IRanges(rep(changepoints$mb_bin_start[i]:changepoints$mb_bin_end[i]))
      overlaps <- IRanges::overlapsAny(kat_range, cp_range)
      if (TRUE %in% overlaps) {
        distance <- 0
      }
      else {
        start_end <- abs(changepoints$mb_bin_start[i] - sample_kat$mb_bin_end[closest_kat])*1e6
        end_start <- abs(changepoints$mb_bin_end[i] - sample_kat$mb_bin_start[closest_kat])*1e6
        distance <- base::min(c(start_end, end_start))
      }
      distances <- c(distances, distance)
    }
    else {
      distances <- c(distances, NA)
    }
  }
  changepoints$distance <- distances
  return (changepoints)
}

cll_changepoints <- findKatDistances('lymph-cll', "CLLE-ES")
cll_changepoints$type <- "Lymph-CLL"
melanoma_changepoints <- findKatDistances('melanoma', c("MELA-AU", "SKCM-US"))
melanoma_changepoints$type <- "Melanoma"
scc_changepoints <- findKatDistances('lung-scc', c("LUSC-CN", "LUSC-US", "LUSC-KR"))
scc_changepoints$type <- "Lung-SCC"
eso_changepoints <- findKatDistances('eso-adenoca', c("ESAD-UK", "ESCA-CN"))
eso_changepoints$type <- "Eso-AdenoCA"
colorect_changepoints <- findKatDistances('colorect', c("COAD-US", "COCA-CN"))
colorect_changepoints$type <- "Colorect-AdenoCA"
bladder_changepoints <- findKatDistances('bladder', c("BLCA-CN", "BLCA-US"))
bladder_changepoints$type <- "Bladder-TCC"
uterus_changepoints <- findKatDistances('uterus-adenoca', c("UTCA-FR", "UCEC-US"))
uterus_changepoints$type <- "Uterus-AdenoCA"
stomach_changepoints <- findKatDistances('stomach', c("STAD-US"))
stomach_changepoints$type <- "Stomach-AdenoCA"
bnhl_changepoints <- findKatDistances('lymph-bnhl', c("DLBC-US", "NHLY-MX"))
bnhl_changepoints$type <- "Lymph-BNHL"
gbm_changepoints <- findKatDistances('cns-gbm', c("GBM-CN", "GBM-US", "LGG-US"))
gbm_changepoints$type <- "CNS-GBM"
bone_changepoints <- findKatDistances('bone-osteosarc', c("BOCA-FR", "BOCA-UK"))
bone_changepoints$type <- "Bone-Osteosarc"
cervix_changepoints <- findKatDistances('cervix', c("CESC-US"))
cervix_changepoints$type <- "Cervix"
head_changepoints <- findKatDistances('head-scc', c("HNCA-MX", "HNSC-US"))
head_changepoints$type <- "Head-SCC"
panc_changepoints <- findKatDistances('panc-adenoca', c("PAAD-US", "PACA-CA","PACA-CN"))
panc_changepoints$type <- "Panc-AdenoCA"
breast_changepoints <- findKatDistances('breast-adenoca', c("BRCA-CN","BRCA-EU","BRCA-FR",
                                                            "BRCA-KR","BRCA-MX","BRCA-UK",
                                                            "BRCA-US"))
breast_changepoints$type <- "Breast-AdenoCA"
prostate_changepoints <- findKatDistances('prost-adenoca', c("PRAD-CA","PRAD-FR","PRAD-CN",
                                                             "PRAD-UK","PRAD-US","PRCA-FR"))
prostate_changepoints$type <- "Prost-AdenoCA"
lung_adeno_changepoints <- findKatDistances('lung-adenoca', "LUAD-US")
lung_adeno_changepoints$type <- "Lung-AdenoCA"
kidney_rcc_changepoints <- findKatDistances('kidney-rcc', c("KIRC-US", "KIRP-US", "RECA-CN",
                                                            "RECA-EU", "CCSK-US"))
kidney_rcc_changepoints$type <- "Kidney-RCC"
kidney_chrcc_changepoints <- findKatDistances('kidney-chrcc', c("KICH-US"))
kidney_chrcc_changepoints$type <- "Kidney-ChRCC"
thyroid_changepoints <- findKatDistances('thy-adenoca', c("THCA-CN", "THCA-US", "THCA-SA"))
thyroid_changepoints$type <- "Thy-AdenoCA"
kataegis_distances <- data.frame()

all_changepoints <- rbind(melanoma_changepoints[,c("distance", "type")], lung_adeno_changepoints[,c("distance", "type")], scc_changepoints[,c("distance", "type")],
                          eso_changepoints[,c("distance", "type")], bladder_changepoints[,c("distance", "type")], colorect_changepoints[,c("distance", "type")],
                          uterus_changepoints[,c("distance", "type")], gbm_changepoints[,c("distance", "type")], cll_changepoints[,c("distance", "type")],
                          bnhl_changepoints[,c("distance", "type")],
                          stomach_changepoints[,c("distance", "type")], breast_changepoints[,c("distance", "type")], panc_changepoints[,c("distance", "type")],
                          kidney_chrcc_changepoints[,c("distance", "type")], kidney_rcc_changepoints[,c("distance", "type")], thyroid_changepoints[,c("distance", "type")],
                          bone_changepoints[,c("distance", "type")], cervix_changepoints[,c("distance", "type")], head_changepoints[,c("distance", "type")],
                          prostate_changepoints[,c("distance", "type")])

peaks <- prostate_peaks
changepoints <- prostate_changepoints
n <- 145
type <- "Prost-AdenoCA"


peaks <- peaks %>%
  dplyr::mutate('end' = start+5,
                'sample_count' = 0)
overlaps <- buildOverlapsData(changepoints, 4)
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
full <- changepoints %>%
  dplyr::mutate(peak = overlaps$peak)
full <- full[,c(1:14,(ncol(full)-2):(ncol(full)))]
mod <- full %>%
  dplyr::filter(peak != 0) %>%
  dplyr::group_by(peak) %>%
  dplyr::summarize(sample_support = dplyr::n_distinct(sampleID)) %>%
  dplyr::filter(sample_support > 2)
full <- full %>%
  tidyr::drop_na() %>%
  dplyr::filter(peak %in% mod$peak) %>%
  dplyr::mutate(type = type)
full <- full[,c(1:14,(ncol(full)-5):ncol(full))]
kataegis_distances <- rbind(kataegis_distances, full)

kataegis_distances$type <- factor(kataegis_distances$type, levels=c("Bladder-TCC",
                                                                    "Lung-SCC",
                                                                    "Melanoma",
                                                                    "Lymph-BNHL",
                                                                    "Bone-Osteosarc",
                                                                    "Cervix",
                                                                    "Eso-AdenoCA",
                                                                    "CNS-GBM",
                                                                    "Stomach-AdenoCA",
                                                                    "Colorect-AdenoCA",
                                                                    "Lymph-CLL",
                                                                    "Uterus-AdenoCA"))

kat_distances_filtered <- kataegis_distances %>%
  dplyr::filter(type != c("Panc-AdenoCA", "Breast-AdenoCA")) %>%
  dplyr::select(-dominant_pre.1, -dominant_post.1)

summary <- kat_distances_filtered %>%
  dplyr::group_by(type, peak) %>%
  dplyr::summarize(median_dist = median(distance),
                   mean_dist = mean(distance),
                   max_dist = max(distance)) %>%
  dplyr::arrange(max_dist)
summary(all_changepoints$distance)
nrow(kat_distances_filtered %>%
       dplyr::filter(distance > 1e6))

ggplot2::ggplot(data=kataegis_distances, ggplot2::aes(x = peak, y = distance, group=peak)) +
  ggplot2::facet_wrap(.~type, nrow=12) +
  ggplot2::geom_vline(data = plyr::ddply(kataegis_distances, "type",
                                         plyr::summarize, peaks = unique(peak)),
                      ggplot2::aes(xintercept=peaks), color='lightblue', alpha=0.5) +
  ggplot2::geom_vline(data = plyr::ddply(kataegis_distances, "type",
                                         plyr::summarize, peaks = unique(peak)+1),
                      ggplot2::aes(xintercept=peaks), color='lightblue', alpha=0.5) +
  ggplot2::geom_vline(data = plyr::ddply(kataegis_distances, "type",
                                         plyr::summarize, peaks = unique(peak)+2),
                      ggplot2::aes(xintercept=peaks), color='lightblue', alpha=0.5) +
  ggplot2::geom_vline(data = plyr::ddply(kataegis_distances, "type",
                                         plyr::summarize, peaks = unique(peak)+3),
                      ggplot2::aes(xintercept=peaks), color='lightblue', alpha=0.5) +
  ggplot2::geom_vline(data = plyr::ddply(kataegis_distances, "type",
                                         plyr::summarize, peaks = unique(peak)+4),
                      ggplot2::aes(xintercept=peaks), color='lightblue', alpha=0.5) +
  #ggplot2::geom_boxplot() +
  ggplot2::geom_jitter(width=2, height=0.1, size=1) +
  ggplot2::theme_bw() +
  ggplot2::scale_x_continuous(breaks = change_bins, labels = chr_labels) +
  ggplot2::scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                              breaks=c(0,10,100,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10)) +
  ggplot2::labs(x = "Chromosome", y = "Distance to Nearest Kataegis Event (bp)") +
  ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 strip.background = ggplot2::element_rect(fill="white"))
ggplot2::ggsave(filename="kataegis_distances.png", plot=ggplot2::last_plot(), device="png", path="~/Desktop/CBSP2021/plots",
                scale=1, width = 30, height = 40, units = "cm", dpi=300)
