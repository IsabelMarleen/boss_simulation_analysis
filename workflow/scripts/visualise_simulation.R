library(ggplot2)
library(argparse)
library(dplyr)
library(patchwork)
library(stringr)

theme_set(theme_minimal())


get_arguments <- function() {
  parser <- argparse::ArgumentParser()
  parser$add_argument('--input_cov', required = TRUE)
  parser$add_argument('--input_unb', required = TRUE)
  parser$add_argument('--dump_time', required = TRUE)
  parser$add_argument('--pores', required = TRUE)
  # parser$add_argument('--genome_size', required = TRUE)
  parser$add_argument('--analysed_log', required = TRUE)
  parser$add_argument('--output', required = TRUE)
  parser$add_argument('--output_seqt', required = TRUE)
  args<- parser$parse_args(commandArgs(trailingOnly = TRUE))
  return(args)
}


visualise_simulation <- function(ptol, nrow) {
  # ARGUMENTS ----
  # get the arguments for input and output files
  args <- get_arguments()
  input_cov <- args$input_cov
  input_unb <- args$input_unb
  dump_time <- as.numeric(args$dump_time)
  pores <- as.numeric(args$pores)
  # genome_size <- as.numeric(args$genome_size)
  log_path <- args$analysed_log
  output <- args$output
  output_seqt <- args$output_seqt
  seq_speed  <- 400

  # ANALYSED LOG ----
  # load the analysed log
  log_file  <- read.csv(log_path) %>%print()

  # extract times for dumps
  log_time  <- log_file %>%
    group_by(cond, dump) %>%
    summarise(time=last(time)) %>%
    rename(seq_time = time, time = dump) %>%print()
 
  # COVERAGE ----
  # load the coverage data
  cov <- read.csv(input_cov) %>%
    left_join(log_time ) %>% 
    mutate(seq_time = if_else(is.na(seq_time), 0, seq_time))

  # cond, time, otu, mean_coverage, low_coverage_prop
  # get the ordering of the OTUs
  otu_order <- cov %>%
    filter(cond == 'control', time == max(time)) %>%
    arrange(desc(mean_coverage)) %>%
    pull(otu)
  # reorder factor levels
  cov$otu <- factor(cov$otu, levels=otu_order)

  # convert time units to minutes
  cov <- cov %>%
  mutate(seq_time = seq_time/seq_speed/pores/60)

  # UNBLOCK ----
  # load the unblocking data and reorder the factor levels
  unb <- read.csv(input_unb) %>%
    left_join(log_time) %>% 
    mutate(orig_time = seq_time) %>%
    mutate(seq_time = if_else(is.na(seq_time), 0, seq_time))
  unb$otu <- factor(unb$otu, levels=otu_order)

  # convert time units to minutes
  unb <- unb %>%
    mutate(seq_time = seq_time/seq_speed/pores/60)

  unb_cond <- unb %>% 
    filter(cond == "boss")


  # convert unblocking data so it is not cumulative
  unb_noncum <- unb %>%
    group_by(cond) %>% #print(n=44)
    mutate(total = total - lag(total, default = first(total))) %>%
    mutate(unb = unb - lag(unb, default = first(unb))) %>%
    mutate(seq_time_noncum = seq_time - lag(seq_time, default = first(seq_time))) %>%
    mutate(base_total = base_total - lag(base_total, default = first(base_total))) %>%
    mutate(unb_ratio = if_else(!is.na(unb/total), unb/total, 0)) %>%
    mutate(total_rate = total/seq_time_noncum)   #%>% #print(n=44)
    # summarise(seq_time=max(seq_time), total=sum(total), unb_ratio=sum(unb)/total, accept_ratio=1-unb_ratio)%>%
    # select(cond,seq_time, unb_ratio, accept_ratio,total) %>% 
    # print(n=44)


  # SANITY CHECKS ----
  # # Calculate coverage manually from all bases sequenced, not from mapping pileup, as a means of checking results
  # unb %>%
  #   mutate(manual_cov = base_total/genome_size)%>%
  #   select(cond, time, manual_cov)%>%
  #   left_join(., select(cov, cond, time, mean_coverage))

  # Quantify time
  alpha = 300
  rho = 300
  mu = 400

  unb_noncum %>%
    mutate(alphas = total * alpha, rhos = unb * rho) %>%
    group_by(cond, time) %>%
    summarise(alphas = sum(alphas), rhos = sum(rhos), time_spent_seq = sum(base_total), dump=max(orig_time, na.rm=T)) %>%
    mutate(total_time = cumsum(alphas) + cumsum(rhos) + cumsum(time_spent_seq), .before = dump) %>%
    mutate(difference = dump - total_time) %>%
    print(n=36)

  unb_noncum %>%
    mutate(alphas = total * alpha, rhos = unb * rho) %>%
    group_by(cond) %>%
    summarise(alphas = sum(alphas), rhos = sum(rhos), time_spent_seq = sum(base_total), dump=max(orig_time, na.rm=T)) %>%
    mutate(total_time = alphas + rhos + time_spent_seq, .before = dump) %>%
    mutate(difference = dump - total_time) %>%
    print()
# 630,919,767
# 565,791,967

# 560056991
# 575991967

  # Coverage dotplots ----
  # load data
  # filenames_boss <- list.files("temp", pattern="*.boss*.csv", full.names=TRUE)
  # ldf_boss <- lapply(filenames_boss, read.csv)
  # cov_boss <- Reduce(full_join, ldf_boss) %>%
  #   pivot_longer(
  #   !V1,
  #   cols_vary = "fastest",
  #   names_to = "t",
  #   values_to = "cov",
  # )

  # filenames_control <- list.files("temp", pattern="*.control*.csv", full.names=TRUE)
  # ldf_control <- lapply(filenames_control, read.csv)
  # cov_control <- Reduce(full_join, ldf_control) %>%
  #   pivot_longer(
  #   !V1,
  #   cols_vary = "fastest",
  #   names_to = "t",
  #   values_to = "cov",
  # )

  # full_cov  <- bind_rows(cov_boss, cov_control)



  # PLOTS ----
  # create plots

  unb_plot <- ggplot(
      data=unb_cond,
      mapping=aes(x=time, y=unb_ratio, color=otu, group=otu), alpha=0.8) +
    geom_line(linewidth=1) +
    geom_point() +
    scale_color_manual(values=ptol) +
    ylab("cum. rejection rate") +
    ylim(0,1)+
    xlim(0, max(unb$seq_time))+
    theme(legend.position = "none")
  # unb_plot

  unb_noncum_plot <- ggplot(
      data=filter(unb_noncum,cond == "boss"),
      mapping=aes(x=time, y=unb_ratio, color=otu, group=otu)) +
    geom_line(linewidth=1) +
    geom_point() +
    scale_color_manual(values=ptol) +
    ylab("rejection rate") +
    ylim(0,1)+
    xlim(0, max(unb$seq_time))+
    theme(legend.position = "none")
  # unb_plot

  nreads <- ggplot(
      data=unb,
      mapping=aes(x=time, y=total, linetype=cond, colour=otu)) +
    geom_line(linewidth=1) +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("# reads")
  # nreads

  nreads_noncum <- ggplot(
     unb_noncum,
      # data=mutate(unb_noncum, total_rate = if_else(time == 0, NA, total_rate)),
      mapping=aes(x=time, y=total_rate, linetype=cond, colour=otu)) +
    geom_line(linewidth=1) +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("# reads/time unit")
  # nreads

  meanc <- ggplot(
      data=cov,
      mapping=aes(x=time, y=mean_coverage, linetype=cond, colour=otu)) +
    geom_line(linewidth=1) +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("mean coverage")
  # meanc

  lowc <- ggplot(
      data=cov,
      mapping=aes(x=time, y=low_coverage_prop, linetype=cond, colour=otu)) +
    geom_line(linewidth=1) +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("prop. sites at <5x")
  # lowc

  evn <- ggplot(
      data=cov,
      mapping=aes(x=time, y=evenness, linetype=cond, colour=otu)) +
    geom_line(linewidth=1) +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("evenness")
  # evenness of coverage

  time_frag <- ggplot(
      data=unb,
      mapping=aes(x=total, y=time, linetype=cond, colour=otu)) +
    geom_line(linewidth=1) +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("seq. time (intervals)")
  # nreads

  # evn_noncum <- ggplot(
  #     data=mutate(cov_noncum),
  #     mapping=aes(x=time, y=evenness, linetype=cond, colour=otu)) +
  #   geom_line(linewidth=1) +
  #   facet_wrap(~otu, scales="free_y", nrow = nrow) +
  #   scale_color_manual(values=ptol, guide = "none") +
  #   ylab("evenness")
  # evenness of coverage

  # create plots with new attempt at scaling time units to real sequencing time

  unb_plot_seq <- ggplot(
      data=unb_cond,
      mapping=aes(x=seq_time, y=unb_ratio, color=otu, group=otu)) +
    geom_line(linewidth=1) +
    geom_point() +
    scale_color_manual(values=ptol) +
    ylab("cum. rejection rate") +
    xlim(0, max(unb$seq_time))+
    ylim(0,1)+
    theme(legend.position = "none")
  # unb_plot

  unb_noncum_plot_seq <- ggplot(
      data=filter(unb_noncum,cond == "boss"),
      mapping=aes(x=seq_time, y=unb_ratio, color=otu, group=otu)) +
    geom_line(linewidth=1) +
    geom_point() +
    ylim(0,1)+
    xlim(0, max(unb$seq_time))+
    scale_color_manual(values=ptol) +
    ylab("rejection rate") +
    theme(legend.position = "none")
  # unb_plot

  nreads_seq <- ggplot(
      data=unb,
      mapping=aes(x=seq_time, y=total, linetype=cond, colour=otu)) +
    geom_line(linewidth=1) +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("# reads")
  # nreads

  nreads_noncum_seq <- ggplot(
      data=mutate(unb_noncum, total_rate = if_else(time == 0, NA, total_rate)),
      mapping=aes(x=seq_time, y=total_rate, linetype=cond, colour=otu)) +
    geom_line(linewidth=1) +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("# reads/min")
  # nreads

  meanc_seq <- ggplot(
      data=cov,
      mapping=aes(x=seq_time, y=mean_coverage, linetype=cond, colour=otu)) +
    geom_line(linewidth=1) +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("mean coverage")
  # meanc

  lowc_seq <- ggplot(
      data=cov,
      mapping=aes(x=seq_time, y=low_coverage_prop, linetype=cond, colour=otu)) +
    geom_line(linewidth=1) +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("prop. sites at <5x")
  # lowc

  evn_seq <- ggplot(
      data=cov,
      mapping=aes(x=seq_time, y=evenness, linetype=cond, colour=otu)) +
    geom_line(linewidth=1) +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("evenness")
  # evenness of coverage

  time_frag_seq <- ggplot(
      data=unb,
      mapping=aes(x=total, y=seq_time, linetype=cond, colour=otu)) +
    geom_line(linewidth=1) +
    facet_wrap(~otu, scales="free_y", nrow = nrow) +
    scale_color_manual(values=ptol, guide = "none") +
    ylab("seq. time (minutes)")
  # nreads


  layout <- ({unb_plot + unb_noncum_plot} / {nreads + nreads_noncum} / {meanc + plot_spacer()} / {lowc + time_frag} / {evn + plot_spacer()} / guide_area()) +
    plot_annotation(tag_levels = "A") +
    plot_layout(guides="collect", heights=rbind(c(1, 1, 1, 1, 1, 0.2),c(1, 1, NA, NA, NA, NA))) &
    xlab("seq. time (intervals)") &
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.text.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )

  layout[[4]] <- layout[[4]] + xlab("# reads")

  ggsave(output, layout, w=8, h=11)

  layout_seqt <- ({unb_plot_seq + unb_noncum_plot_seq} / {nreads_seq + nreads_noncum_seq} / {meanc_seq + time_frag_seq} / {lowc_seq + plot_spacer()} / {evn_seq + plot_spacer()}/ guide_area()) +
    plot_annotation(tag_levels = "A") +
    plot_layout(guides="collect", heights=rbind(c(1, 1, 1, 1, 1, 0.2),c(1, 1, 1, NA, NA, NA))) &
    xlab("seq. time (minutes)") &
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.text.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # layout_seqt <- wrap_plots(list(unb_plot_seq, lowc_seq, nreads_seq, plot_spacer(), meanc_seq, plot_spacer()), axes='keep', axis_titles="keep")+
  #   plot_annotation(tag_levels = "A") +
  #   plot_layout(
  #     nrow = 3,
  #     ncol = 2,
  #     guides="collect"
  #     ) &
  #   xlab("seq. time (pseudominutes)") &
  #   theme(
  #     legend.position = "bottom",
  #     legend.title = element_blank(),
  #     strip.text.x = element_blank(),
  #     plot.background = element_rect(fill = "white", color = NA)
  #   )

  # layout_seqt[[6]] <- layout_seqt[[3]] + xlab("# reads")
  ggsave(output_seqt, layout_seqt, w=8, h=7)

}


ptol <- c("#CC6677","#332288","#DDCC77","#117733","#88CCEE","#882255","#44AA99","#999933","#AA4499")
nrow <- 1
visualise_simulation(ptol, nrow)



