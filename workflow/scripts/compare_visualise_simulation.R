library(ggplot2)
library(argparse)
library(dplyr)
library(patchwork)

theme_set(theme_minimal())


ptol <- c("#CC6677","#332288","#DDCC77","#117733","#88CCEE","#882255","#44AA99","#999933","#AA4499")

nrow <- 1


get_arguments <- function() {
  parser <- argparse::ArgumentParser()
  parser$add_argument('--input_cov_A', required = TRUE)
  parser$add_argument('--input_unb_A', required = TRUE)
parser$add_argument('--input_cov_B', required = TRUE)
  parser$add_argument('--input_unb_B', required = TRUE)
  parser$add_argument('--dump_time', required = TRUE)
  parser$add_argument('--pores', required = TRUE)
  parser$add_argument('--genome_size', required = TRUE)
  parser$add_argument('--analysed_log', required = TRUE)
  parser$add_argument('--analysed_log2', required = TRUE)
  parser$add_argument('--output', required = TRUE)
  parser$add_argument('--output_seqt', required = TRUE)
  args<- parser$parse_args(commandArgs(trailingOnly = TRUE))
  return(args)
}


# get the arguments for input and output files
args <- get_arguments()
input_cov_a <- args$input_cov_A
input_unb_a <- args$input_unb_A
input_cov_b <- args$input_cov_B
input_unb_b <- args$input_unb_B
dump_time <- as.numeric(args$dump_time)
pores <- as.numeric(args$pores)
genome_size <- as.numeric(args$genome_size)
log_a <- args$analysed_log
log_b <- args$analysed_log2
output <- args$output
output_seqt <- args$output_seqt

dump_time1 <- 35000000
seq_speed  <- 400

# ANALYSED LOG ----
# load the analysed log
log_file_a  <- read.csv(log_a) %>%print()

# extract times for dumps
log_time_a  <- log_file_a %>%
  group_by(cond, dump) %>%
  summarise(time=last(time)) %>%
  rename(seq_time = time, time = dump) %>%print()
 
# load the analysed log
log_file_b  <- read.csv(log_b) %>%print()

# extract times for dumps
log_time_b  <- log_file_b %>%
  group_by(cond, dump) %>%
  summarise(time=last(time)) %>%
  rename(seq_time = time, time = dump) %>%print()

# load the coverage data
cov_a <- read.csv(input_cov_a) %>%
  left_join(log_time_a ) %>% 
  mutate(seq_time = if_else(is.na(seq_time), 0, seq_time))

cov_b <- read.csv(input_cov_b) %>%
    left_join(log_time_b ) %>% 
    mutate(seq_time = if_else(is.na(seq_time), 0, seq_time))

cov <- bind_rows('accept unmapped'=cov_a, 'reject unmapped'=cov_b, .id="sim")
# cond, time, otu, mean_coverage, low_coverage_prop
# get the ordering of the OTUs
otu_order <- cov %>%
  filter(cond == 'control', time == max(time)) %>%
  arrange(desc(mean_coverage)) %>%
  distinct(otu) %>%
  pull(otu)
# reorder factor levels
cov$otu <- factor(cov$otu, levels=otu_order)


# load the unblocking data and reorder the factor levels
unb <- bind_rows('accept unmapped'=read.csv(input_unb_a), 'reject unmapped'=read.csv(input_unb_b), .id="sim")%>%
  left_join(bind_rows('accept unmapped'=log_time_a, 'reject unmapped'=log_time_b, .id="sim")) %>%
  mutate(seq_time = if_else(is.na(seq_time), 0, seq_time))
unb$otu <- factor(unb$otu, levels=otu_order)

unb <- unb %>%
  mutate(seq_time = seq_time/pores/60/seq_speed) # convert time units to minutes

unb_cond <- unb %>% 
  filter(cond == "boss")


cov <- cov %>%
  mutate(seq_time = seq_time/pores/60/seq_speed)  # convert time units to minutes

# convert unblocking data so it is not cumulative
unb_noncum <- unb %>%
  group_by(cond, sim) %>%
  mutate(total = total - lag(total, default = first(total))) %>%
  mutate(unb = unb - lag(unb, default = first(unb))) %>%
  mutate(base_total = base_total - lag(base_total, default = first(base_total))) %>%
  mutate(unb_ratio = if_else(!is.na(unb/total), unb/total, 0))

# create plots with new attempt at scaling time units to real sequencing time

unb_plot_seq <- ggplot(
    data=unb_cond,
    mapping=aes(x=seq_time, y=unb_ratio, color=sim, group=sim)) +
  geom_line(linewidth=1) +
  geom_point() +
  scale_color_manual(values=ptol) +
  ylab("cum. rejection rate") +
  ylim(0,1)+
  xlim(0, max(unb$seq_time))+
  theme(legend.position = "none")
# unb_plot

unb_noncum_plot_seq <- ggplot(
    data=filter(unb_noncum,cond == "boss"),
    mapping=aes(x=seq_time, y=unb_ratio, color=sim)) +
  geom_line(linewidth=1) +
  geom_point() +
  scale_color_manual(values=ptol) +
  ylab("rejection rate") +
  ylim(0,1)+
  xlim(0, max(unb$seq_time))+
  theme(legend.position = "none")
# unb_plot

nreads_seq <- ggplot(
    data=unb,
    mapping=aes(x=seq_time, y=total, linetype=cond, colour=sim)) +
  geom_line(linewidth=1) +
  facet_wrap(~otu, scales="free_y", nrow = nrow) +
  scale_color_manual(values=ptol, guide = "none") +
  ylab("# reads")
# nreads

nreads_noncum_seq <- ggplot(
    data=mutate(unb_noncum, total = if_else(time == 0, NA, total)),
    mapping=aes(x=seq_time, y=total, linetype=cond, colour=sim)) +
  geom_line(linewidth=1) +
  facet_wrap(~otu, scales="free_y", nrow = nrow) +
  scale_color_manual(values=ptol, guide = "none") +
  ylab("# reads")
# nreads

meanc_seq <- ggplot(
    data=cov,
    mapping=aes(x=seq_time, y=mean_coverage, linetype=cond, colour=sim)) +
  geom_line(linewidth=1) +
  facet_wrap(~otu, scales="free_y", nrow = nrow) +
  scale_color_manual(values=ptol, guide = "none") +
  ylab("mean coverage")
# meanc

lowc_seq <- ggplot(
    data=cov,
    mapping=aes(x=seq_time, y=low_coverage_prop, linetype=cond, colour=sim)) +
  geom_line(linewidth=1) +
  facet_wrap(~otu, scales="free_y", nrow = nrow) +
  scale_color_manual(values=ptol, guide = "none") +
  ylab("prop. sites at <5x")
# lowc



layout_seqt <- wrap_plots(list(unb_plot_seq, unb_noncum_plot_seq, nreads_seq, nreads_noncum_seq, meanc_seq, lowc_seq), axes='keep', axis_titles="keep")+
  plot_annotation(tag_levels = "A") +
  plot_layout(
    nrow = 3,
    ncol = 2,
    guides="collect"
    ) &
  xlab("seq. time (pseudominutes)") &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text.x = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

# layout_seqt[[6]] <- layout_seqt[[3]] + xlab("# reads")
ggsave(output_seqt, layout_seqt, w=8, h=7)