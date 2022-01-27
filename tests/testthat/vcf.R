cll_vcfs <- list()
sbs1 <- c()
sbs2 <- c()
sbs5 <- c()
sbs7a <- c()
sbs7b <- c()
sbs7c <- c()
sbs7d <- c()
sbs13 <- c()
sbs17a <- c()
sbs17b <- c()
sbs38 <- c()
sbs40 <- c()
sbs9 <- c()
for (i in c(1:330)) {
  if (file.exists(paste0("~/Desktop/CBSP2021/",'cll_vcf/',"cll_vcf",i,".Rdata"))) {
    print(i)

    load(paste0("~/Desktop/CBSP2021/",'cll_vcf/',"cll_vcf",i,".Rdata"))
    temp <- load(paste0("~/Desktop/CBSP2021/",'cll_vcf/',"cll_vcf",i,".Rdata"))
    sbs1 <- c(sbs1, ((get(temp)[['mixtures']][1,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][1,1])))
    #sbs2 <- c(sbs2, ((get(temp)[['mixtures']][2,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][2,1])))
    sbs5 <- c(sbs5, ((get(temp)[['mixtures']][2,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][2,1])))
    #sbs7a <- c(sbs7a, ((get(temp)[['mixtures']][4,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][4,1])))
    #sbs7b <- c(sbs7b, ((get(temp)[['mixtures']][5,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][5,1])))
   # sbs7c <- c(sbs7c, ((get(temp)[['mixtures']][6,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][6,1])))
   # sbs7d <- c(sbs7d, ((get(temp)[['mixtures']][7,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][7,1])))
    sbs9 <- c(sbs9, ((get(temp)[['mixtures']][3,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][3,1])))

  #  sbs13 <- c(sbs13, ((get(temp)[['mixtures']][8,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][8,1])))
   # sbs17a <- c(sbs17a, ((get(temp)[['mixtures']][9,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][9,1])))
  #  sbs17b <- c(sbs17b, ((get(temp)[['mixtures']][10,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][10,1])))
  #  sbs38 <- c(sbs38, ((get(temp)[['mixtures']][11,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][11,1])))
    sbs40 <- c(sbs40, ((get(temp)[['mixtures']][4,ncol(get(temp)[['mixtures']])-1])-(get(temp)[['mixtures']][4,1])))

  }
}
cll_palette <- default_palette[c(1,5,10,21)]
cll_changes <- data.frame("SBS1" = sbs1, "SBS5" = sbs5, "SBS9" = sbs9, "SBS40" = sbs40)
                               # "SBS7a" = sbs7a, "SBS7b" = sbs7b, "SBS7c" = sbs7c,
                               # "SBS7d" = sbs7d, "SBS13" = sbs13, "SBS17a" = sbs17a,
                               # "SBS17b" = sbs17b, "SBS38" = sbs38, "SBS40"=sbs40)
cll_changes <- cll_changes %>%
  tidyr::pivot_longer(cols=c("SBS1", "SBS5", "SBS9", "SBS40"), names_to="Signature", values_to="exposure_change")
cll_changes$Signature <- factor(cll_changes$Signature, levels=c("SBS1", "SBS5", "SBS9", "SBS40"))
ggplot2::ggplot(data=cll_changes, ggplot2::aes(x=Signature, y = exposure_change*100, fill=Signature)) +
  ggplot2::geom_violin() +
  ggplot2::scale_fill_manual(values=cll_palette) +
  ggplot2::scale_color_manual(values=cll_palette) +
  ggplot2::scale_y_continuous(limits=c(-80,80), breaks=c(-80,-40,0,40,80), labels=c(-80,-40,0,40,80)) +
  ggplot2::theme_bw() +
  ggplot2::labs(y = "Exposure Change (%)") +
 ggplot2::stat_summary(fun.data=ggplot2::mean_sdl, mult=1,
                  geom="pointrange", color="black") +
  ggplot2::theme(legend.position = 'none')

ggplot2::ggsave(filename="cll_violins.pdf", path="~/Desktop/CBSP2021/plots", plot=ggplot2::last_plot(), device="pdf",
                scale=1, width=20, height=15, unit="cm")


p2 <- plotSpaceTrajectory(traj, chr_level=T, cutoff=0) +
  ggplot2::scale_color_manual(values=lung_scc_palette)
p1 <- plotSpaceTrajectory(lung6, chr_level=F, cutoff=0) +
  ggplot2::scale_color_manual(values=lung_scc_palette) +
  ggplot2::scale_y_continuous(breaks=c(0,20,40,60,80), labels=c(0,20,40,60,80), limits=c(0,80))


prow <- cowplot::plot_grid(
  p1 + ggplot2::theme(legend.position="none"),
  p2 + ggplot2::theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)

# extract the legend from one of the plots
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  p1 + ggplot2::theme(legend.box.margin = ggplot2::margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of
# the width of one plot (via rel_widths).
prow <- cowplot::plot_grid(prow, legend, rel_widths = c(6, .4))

full_plot <- cowplot::plot_grid(p3, prow, labels = c('A', ''), label_size = 12, ncol = 1)

ggplot2::ggsave(filename="lung-scc6_pairwise.png", path="~/Desktop/CBSP2021/plots", plot=ggplot2::last_plot(), device="png",
                scale=1, width=50, height=10, unit="cm")


ggplot2::ggplot(data = changepoint_positions, ggplot2::aes(x=changepoint_positions, color="red")) +
  ggplot2::geom_histogram(fill="red", binwidth=2)
  ggplot2::theme_bw() +
  ggplot2::scale_y_continuous(limits=c(0,0.1))
  ggplot2::labs(x = "Chromosome", y = "Changepoint Frequency", fill = 'Distribution') +
  ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank())

melanoma_palette <- c(default_palette[1],default_palette[2],default_palette[5],
                      default_palette[7],default_palette[8],"orange",
                      "lightgreen",default_palette[14],default_palette[17],
                      default_palette[18],"salmon",default_palette[21])


for (i in c(1:107)) {
  if (file.exists(paste0("~/Desktop/CBSP2021/melanoma_sigfreq/melanoma_sigfreq",i,".Rdata"))) {
    load(paste0("~/Desktop/CBSP2021/melanoma_sigfreq/melanoma_sigfreq",i,".Rdata"))
    temp <- load(paste0("~/Desktop/CBSP2021/melanoma_sigfreq/melanoma_sigfreq",i,".Rdata"))

    binData <- traj[['binData']] %>%
      dplyr::mutate(start_chrom=chr,
                    end_chrom = chr)
    binData <- binData %>%
      dplyr::select(start_chrom, end_chrom, pos, bin, VAF) %>%
      dplyr::group_by(bin) %>%
      dplyr::summarize(start_chrom = min(start_chrom),
                       end_chrom = max(end_chrom),
                       start = min(pos),
                       end = max(pos),
                       mean_vaf = mean(VAF))
    for (i in 1:(nrow(binData)-1)) {
      if (binData$end_chrom[i] > binData$start_chrom[i]) {
        binData$end[i] <- binData$start[i+1]-1
      }
    }
    traj[[4]] <- binData
    colnames(traj[['mixtures']]) <- binData$bin
    p1 <- plotSpaceTrajectory(traj, chr_level=F, cutoff=0)
    p2 <- ggplot2::ggplot(data = binData, ggplot2::aes(x = bin, y = mean_vaf)) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::labs(y = "VAF") +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank())

    plots <- cowplot::align_plots(p1, p2, align='v', axis='tblr')


    prow <- cowplot::plot_grid(
      plots[[2]],
      plots[[1]],
      ncol=1,
      rel_heights = c(0.5,4)
    )
    prow
    ggplot2::ggsave(filename=paste0("melanoma_sigfreq",55,".png"), path="~/Desktop/CBSP2021/plots", plot=ggplot2::last_plot(), device="png",
                    scale=1, width=40, height=20, unit="cm")
  }
}
