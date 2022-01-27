library(reticulate)
np <- import("numpy")

dataframe <- np$load("data_SKCM_w1e6.npz")
hist_counts <- data.table::transpose(as.data.frame(dataframe$f[["hist_counts"]]))
targets <- data.table::transpose(as.data.frame(dataframe$f[["targets"]]))

data <- data.frame("hist_counts" = hist_counts, "targets" = targets)
colnames(data) <- c("hist_counts", 'targets')
data <- data[1:3113,]

calculateDNAseData <- function(master, binSize) {
  binned_master <- getBinNumber(master, binSize)
  data$bin <- binned_master$bin

  dataframe_sum <- data %>%
    dplyr::group_by(bin) %>%
    dplyr::summarize(min_count = min(hist_counts),
                     mean_count = mean(hist_counts),
                     max_count = max(hist_counts))

  return (dataframe_sum)
}

calculateMutDensity <- function(trajectory) {
  width <- bin <- mut_density <- NULL

  binData <- trajectory[[4]]
  # calculate raw mutation density at each bin as # mutations / 1 bp
  widths <- binData %>%
    dplyr::mutate(rowsum = rowSums(binData[7:102])) %>%
    dplyr::select(width, rowsum) %>%
    dplyr::mutate(mut_density = rowsum / width)

  # assign the correct data column to bin number-- depends on SpaceTrack() parameters
  if (!is.null(binData$genome_bin)) {
    widths$bin <- binData$genome_bin
  }
  else if (!is.null(binData$actual_bin)) {
    widths$bin <- binData$actual_bin
  }
  else {
    widths$bin <- binData$bin
  }

  #apply Min-Max normalization to mut_density column
  density_norm <- sapply(widths[3], min_max_norm)

  widths <- widths %>%
    dplyr::mutate(mut_density = density_norm)

  return (widths)
}

# extract gene density data for human genome from karyoPloteR gene density ideogram
txdb <-TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
all.genes <- base::suppressMessages(GenomicFeatures::genes(txdb))
chrom_lengths <- stats::setNames(object = all.genes@seqinfo@seqlengths[1:24], all.genes@seqinfo@seqnames[1:24])
windows <- GenomicRanges::tileGenome(seqlengths = chrom_lengths,
                                     tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
dens <- GenomicRanges::countOverlaps(windows, all.genes)

calculateGeneDensity <- function(master, binSize, trajectory){

  start_chrom <- bin <- density <- gene_density <- mean_density <- NULL

  # assign bins to 1Mb regions of genome
  binned_master <- getBinNumber(master, binSize)

  # build dataframe of gene density values for each bin

  density_data <- data.frame(start_chrom = binned_master$start_chrom,
                             bin = binned_master$bin, density = dens)

  if (!is.null(trajectory[[4]]$genome_bin)) {
    # if plottting a shuffle trajectory, exclude Y chromosome data
    density_data <- density_data %>%
      dplyr::filter(start_chrom != 24) %>%
      dplyr::group_by(bin) %>%
      dplyr::summarize(mean_density = base::mean(density))

  }
  else {
    density_data <- density_data %>%
      dplyr::group_by(bin) %>%
      dplyr::summarize(mean_density = base::mean(density))
  }

  #apply Min-Max normalization to gene density column
  density_norm <- sapply(density_data[2], min_max_norm)

  density_data <- density_data %>%
    dplyr::mutate(mean_density = density_norm)

  return (density_data)
}

file_ids <- read.delim(paste0("~/Desktop/CBSP2021/",'melanoma',"_ids.txt"), header=FALSE)
colnames(file_ids) <- c("id")

corrs <- data.frame("id" = c(), "n_cps" = c(), "dnase_corr" = c(), "mut_corr" = c(), "gene_corr" = c())

for (i in c(1:nrow(file_ids))) {
  fileID <- trimws(as.character(file_ids$id[52]), which='both')
  if (file.exists(paste0("~/Desktop/CBSP2021/melanoma/melanoma", i, ".Rdata"))) {
    master <- readFormat(paste0("~/Desktop/CBSP2021/melanoma_masters/",fileID,".MBcounts.csv"))

    load(paste0("~/Desktop/CBSP2021/",'melanoma',"/","melanoma",i,".Rdata"))
    temp <- load(paste0("~/Desktop/CBSP2021/",'melanoma',"/","melanoma",i,".Rdata"))

    mutDensity <- calculateMutDensity(get(temp))
    if (mutDensity$bin[1] != 1) {
      mutDensity <- mutDensity[order(nrow(mutDensity):1),]
    }
    geneDensity <- calculateGeneDensity(master, 200, get(temp))
    dnasePeaks <- calculateDNAseData(master, 200)
    features_data <- data.frame("id" = i, "bin" = dnasePeaks$bin, "dnase_peaks" = dnasePeaks$mean_count,
                                "gene_density" = geneDensity$mean_density,
                                "mut_density" = mutDensity$mut_density)
    cpPos <- assignChangepoints(get(temp), 0)
    features_data <- features_data %>%
      dplyr::mutate(cp = dplyr::if_else(features_data$bin %in% cpPos$cpPos1,
                                        1, 0))

    dnase_model <- glm(cp ~ dnase_peaks, family=binomial(link='logit'), data=features_data)
    dnase_cor <- sqrt(pscl::pR2(dnase_model)[4])

    mut_model <- glm(cp ~ mut_density, family=binomial(link='logit'), data=features_data)
    mut_cor <- sqrt(pscl::pR2(mut_model)[4])

    gene_model <- glm(cp ~ mean_density, family=binomial(link='logit'), data=features_data)
    gene_cor <- sqrt(pscl::pR2(gene_model)[4])

    corrs <- rbind(corrs, c(features_data$id[1], nrow(cpPos), dnase_cor, mut_cor, gene_cor))
  }
}

colnames(corrs) <- c("id", "n_cps", "dnase_cor", "mut_cor", "gene_cor")
corrs <- corrs %>%
  dplyr::filter(n_cps > 0)


p1 <- plotGenomicFeaturesTrajectory(master, traj, 200, chr_level=T)
ggplot2::ggsave(filename="melanoma52_full.png", path="~/Desktop/CBSP2021/plots", device="png",
       plot=ggplot2::last_plot(), scale=1, width=34, height=24, unit="cm", dpi="retina")


p4 <- ggplot2::ggplot(data=corrs, ggplot2::aes(x=gene_cor)) +
  ggplot2::geom_histogram(binwidth=0.01, fill='cornflowerblue') +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Correlation") +
  ggplot2::scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                              limits=c(0,1)) +
  ggplot2::scale_y_continuous(breaks=c(0,2,4,6), labels=c(0,2,4,6), limits=c(0,7)) +
  ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
  ggplot2::annotate("rect", xmin=mean(corrs$gene_cor)-0.005, xmax=mean(corrs$gene_cor)+0.005,
                    ymin=-Inf, ymax=Inf, fill='red')
p4

prow <- cowplot::plot_grid(
  p2,
  p4,
  p3,
  align = 'vh',
  labels = c("B", "C", "D"),
  label_size=10,
  hjust = -1,
  nrow = 1
)
prow2 <- cowplot::plot_grid(
  p1,
  align = 'vh',
  labels = c("A"),
  label_size=10,
  hjust = -1,
  nrow = 1
)
# extract the legend from one of the plots
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  p1 + ggplot2::theme(legend.box.margin = ggplot2::margin(0, 0, 0, 12))
)
prow3 <- cowplot::plot_grid(prow2, legend, rel_widths = c(6, .4))

full_plot <- cowplot::plot_grid(prow2, prow, label_size = 10, ncol=1, rel_heights = c(1,0.3))
full_plot
ggplot2::ggsave(filename="signature-tracks.pdf", path="~/Desktop/CBSP2021/plots", device="pdf",
                plot=ggplot2::last_plot(), scale=1, width=50, height=25, unit="cm", dpi="retina")
