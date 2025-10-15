library(ggplot2)
library(argparse)
library(dplyr)
library(patchwork)
library(stringr)
library(readr)
library(languageserversetup)
# library(languageserver)

# languageserver_add_to_rprofile()

# languageserver_install()

get_arguments <- function() {
    parser <- argparse::ArgumentParser()
    parser$add_argument('--input_accept', required = TRUE)
    parser$add_argument('--input_reject', required = TRUE)
    parser$add_argument('--output', required = TRUE)
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return(args)
}

analyse_log  <- function() {
    # Obtained by running the following in a terminal:
    # awk -v var=$EXP 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths) {print l, lengths[l]}; print "total reads: " counter}' /nfs/research/goldman/ipoetzsch/trio_exp2/nanosim_out_2mil/simulated_2mil_filt_aligned_reads.fq | sort -nr > /hps/nobackup/goldman/ipoetzsch/length_distr_nanosim_2mil.txt
    
    orig_read_lengths_path  <- "/hps/nobackup/goldman/ipoetzsch/length_distr_nanosim_2mil.txt"
    orig_read_lengths <- 
    read_table(orig_read_lengths_path, comment = "total", col_names=c("length","count")) %>% arrange(length) #%>%print()

    orig_read_lengths_path_no_filter <- "/hps/nobackup/goldman/ipoetzsch/length_distr_nanosim_2mil_combined.txt"
    orig_read_lengths_no_filter <- 
    read_table(orig_read_lengths_path_no_filter, comment = "total", col_names=c("length","count")) %>% arrange(length) #%>%print()

    args <- get_arguments()
    input_accept <- args$input_accept
    input_reject <- args$input_reject
    output <- args$output

    header <- c("read_id", "fate_contr", "length_contr", "map_trunc", "map_full", "fate", "length")
    accept <- read_csv(input_accept, col_names=header)
    reject <- read_csv(input_reject, col_names=header)

    joint <- bind_rows(list(accept = accept, reject = reject), .id='exp')
    joint %>%
        group_by(exp) %>%
        summarise(no_map_trunc = 200000-sum(map_trunc), no_map_full = 200000 - sum(map_full), rejections = 200000-sum(fate), 
        length=sum(length), avg_length = mean(length_contr), length_contr = sum(length_contr), boss_time = length/512/400/60) %>%
        mutate(true_rejections = if_else(exp=="reject", rejections - no_map_trunc, rejections)) %>%
        mutate(diff_in_rejections = Reduce("-", rejections, accumulate = T)) %>%
        mutate(cond_diff_length = Reduce("-", length, accumulate = T))%>%
        mutate(diff_true_rejections = Reduce("-", true_rejections, accumulate = T))%>%
        mutate(length_diff = length_contr - length)%>%
        print(width=Inf)

    # Get reads that did not map and were read in accept but not in reject and get the difference in length
    joint %>%
        group_by(exp) %>%
        filter(map_trunc==0)%>% #print()
        summarise(
            length_contr_sum = sum(length_contr), length_sum = sum(length), avg_contr = mean(length_contr), 
            avg = mean(length), sd_contr = sd(length_contr), sd =sd(length))%>%
        mutate(diff = length_contr_sum - length_sum, .before=avg_contr)%>%
        # length=sum(length), avg_length = mean(length_contr)) %>%
        # mutate(no_map = no_map_trunc+no_map_full)%>%
        # mutate(diff_in_accept = rejections - no_map) %>%
        # mutate(diff_in_rejections = Reduce("-", rejections, accumulate = T)) %>%
        print(width=Inf)
    joint %>%
        group_by(exp)%>%
        filter(length==400)%>%
        summarise(rej = n())%>%print()

    # Plot fragment length distribution for the three conditions
    length_dist <- joint %>%
        ggplot()+
        geom_vline(aes(xintercept=400), linetype="dotted")+
        geom_freqpoly(aes(length_contr, after_stat(count), group=exp, colour="control", linetype="control"), binwidth=200)+
        geom_freqpoly(aes(length, after_stat(count), group=exp, colour=exp, linetype=exp), alpha=0.6,binwidth=200)+
        scale_y_continuous(trans = 'log10')+
        ylab("Fragment count")+
        xlab("Fragment length (bp)")+
        labs(color = "Experiment", linetype = "Experiment")+     
        scale_linetype_manual(name="Experiment", values=c(1,1,2)) #+
        # scale_linetype_manual(name="Experiment", values=c(1,2))
        # scale_color_manual(name="Experiment", values=c(3)) 
        # guides(color="none",linetype="none")

    length_dist_density <- joint %>%
        ggplot()+
        geom_vline(aes(xintercept=400), linetype="dotted")+
        geom_freqpoly(aes(length_contr, after_stat(density), group=exp, colour="control", linetype="control"), binwidth=800)+
        geom_freqpoly(aes(length, after_stat(density), group=exp, colour=exp, linetype=exp), alpha=0.6,binwidth=800)+
        scale_y_continuous(trans = 'log10')+
        ylab("Fragment density")+
        xlab("Fragment length (bp)")+
        # ylim(0,1)+
        labs(color = "Experiment", linetype = "Experiment")+     
        scale_linetype_manual(name="Experiment", values=c(1,1,2)) #+
        # scale_linetype_manual(name="Experiment", values=c(1,2))
        # scale_color_manual(name="Experiment", values=c(3)) 
        # guides(color="none",linetype="none")
    ggsave(output, length_dist_density)
    
    length_dist_density_input  <- orig_read_lengths%>%
        ggplot()+
        geom_vline(aes(xintercept=400), linetype="dotted")+
        geom_line(aes(length, count, color="> 401bp", linetype="> 401bp"), stat="summary_bin", binwidth = 200, fun=sum)+
        geom_line(data=orig_read_lengths_no_filter, aes(length, count, color="no filter", linetype="no filter"), stat="summary_bin", binwidth = 200, fun=sum)+
        geom_freqpoly(data=joint, aes(length_contr, after_stat(count), group=exp, colour="> 401bp, subset", linetype="> 401bp, subset"), binwidth=200)+
        # geom_line(aes(length, count))+
        # geom_freqpoly(aes(length, after_stat(density), group=exp, colour=exp, linetype=exp), alpha=0.6,binwidth=800)+
        # scale_y_continuous(trans="log10")+
        ylab("Fragment count")+
        xlab("Fragment length (bp)")+
        # ylim(0,1)+
        # labs(color = "Experiment", linetype = "Experiment")+     
        scale_linetype_manual(name="Experiment", values=c(1,1,2)) +
        # scale_linetype_manual(name="Experiment", values=c(1,2))+
        scale_color_discrete(name="Experiment", values=c(1,2,3)) 
        # guides(color="none",linetype="none")
    ggsave(output, length_dist_density_input)

    length_dist_density_input_no_filter  <- orig_read_lengths_no_filter%>%
        ggplot()+
        geom_vline(aes(xintercept=400), linetype="dotted")+
        geom_line(aes(length, count), stat="summary_bin", binwidth = 200, fun=sum)+
        # geom_line(aes(length, count))+
        # geom_freqpoly(aes(length, after_stat(density), group=exp, colour=exp, linetype=exp), alpha=0.6,binwidth=800)+
        # scale_y_continuous(trans = 'log10')+
        ylab("Fragment count")+
        xlab("Fragment length (bp)")
        # ylim(0,1)+
        # labs(color = "Experiment", linetype = "Experiment")+     
        # scale_linetype_manual(name="Experiment", values=c(1,1,2)) #+
        # scale_linetype_manual(name="Experiment", values=c(1,2))
        # scale_color_manual(name="Experiment", values=c(3)) 
        # guides(color="none",linetype="none")
    ggsave(output, length_dist_density_input_no_filter)

    length_dist <- joint %>%
        ggplot()+
        geom_vline(aes(xintercept=400), linetype="dotted")+
        geom_freqpoly(aes(length_contr, after_stat(count), group=exp, colour="control", linetype="control"), binwidth=1000)+
        geom_freqpoly(aes(length, after_stat(count), group=exp, colour=exp, linetype=exp), alpha=0.6,binwidth=1000)+
        # scale_y_continuous(trans = 'log10')+
        ylab("Fragment count")+
        xlab("Fragment length (bp)")+
        # ylim(0,1)+
        ylim(0,200000)+
        labs(color = "Experiment", linetype = "Experiment")+     
        scale_linetype_manual(name="Experiment", values=c(1,1,2)) #+
        # scale_linetype_manual(name="Experiment", values=c(1,2))
        # scale_color_manual(name="Experiment", values=c(3)) 
        # guides(color="none",linetype="none")
    ggsave(output, length_dist)

    length_dist_contr <- joint %>%
        ggplot()+
        geom_vline(aes(xintercept=400), linetype="dotted")+
        geom_freqpoly(aes(length_contr, after_stat(count), group=exp, colour="control", linetype="control"), binwidth=1)+
        # geom_freqpoly(aes(length, after_stat(count), group=exp, colour=exp, linetype=exp), alpha=0.6,binwidth=1)+
        # scale_y_continuous(trans = 'log10')+
        ylab("Fragment count")+
        xlab("Fragment length (bp)")+
        # ylim(0,1)+
        labs(color = "Experiment", linetype = "Experiment")+     
        # scale_linetype_manual(name="Experiment", values=c(1,1,2)) #+
        scale_linetype_manual(name="Experiment", values=c(1,2))
        # scale_color_manual(name="Experiment", values=c(3)) 
        # guides(color="none",linetype="none")

     ggsave(output, length_dist)

    d_contr <- ggplot_build(length_dist_contr)
    d <- ggplot_build(length_dist)
d$data[[2]] %>% 
        bind_rows(d_contr$data[[2]], .id="exp") %>%
        mutate(exp = if_else(exp==1, "exp", "control"))%>%
        filter(!(exp=="control" & group==2))%>%
        mutate(exp = if_else(exp=="exp"&group==1, "accept", exp))%>%
        mutate(exp = if_else(group==2, "reject", exp)) %>%
        select(x,y, count, ncount, exp) %>% #head()%>%print()
        mutate(len_group = if_else(x>400, ">400", "400")) %>%#distinct(exp, len_group)%>%print()
        group_by(len_group, exp) %>%
        dplyr::summarise(y_sum=sum(y), count_sum=sum(count), ncount_sum=sum(ncount), .groups="keep") %>% #print()
        # mutate(total_reads = if_else(exp))
        mutate(ratio = count_sum/200000)%>%
        print()

    d_contr$data[[2]] %>%
        filter(group == 1) %>% # contains data twice for each exp group, but they are identical, see table(select(filter(d_contr$data[[1]], group==1), -group) == select(filter(d_contr$data[[1]], group==2),-group))
        select(x,y, count, ncount) %>%
        # mutate(total=sum(y)) %>%print()
        mutate(len_group = if_else(x>400, ">400", "400")) %>%
        mutate(total = sum(y))%>%
        group_by(len_group) %>%
        summarise(y=sum(y), count=sum(count), ncount=sum(ncount)) %>%
        # mutate(total_reads = if_else(exp))
        mutate(ratio = count/200000)%>%
        print()

    # Length distributions of reads that don't map
    joint %>%
        group_by(map_full, map_trunc, exp)%>%
        summarise(avg_length_contr = mean(length_contr), sd_length_contr = sd(length_contr), count=n(), avg_length = mean(length), sd_length = sd(length))%>%
        # select(count, avg_length, sd_length)%>%
        print(width=Inf)    

    # Plot fragment length distribution for the reads that don't map
    length_dist_nomap <- joint %>%
        filter(map_trunc==0)%>%
        ggplot()+
        geom_freqpoly(aes(length_contr, after_stat(count), group=exp, colour=exp, linetype=exp),binwidth=800)+
        geom_vline(aes(xintercept=400), linetype="dotted")+
        ylab("Fragment count")+
        xlab("Fragment length (bp)")+
        labs(color = "Experiment", linetype = "Experiment")+     
        scale_linetype_manual(values=c(1,2))
    ggsave("/hps/nobackup/goldman/ipoetzsch/boss_simulation_example/results/comp_accept_reject/2025-07-17_log_reads_accept_reject_comp_dist_nomap.png", length_dist_nomap)    
    
    # Plot fragment length distribution split by reads that do or don't map at 400bp
    length_dist_by_map <- joint %>%
        ggplot()+
        geom_freqpoly(aes(length_contr, after_stat(count), group=interaction(map_trunc,exp), colour=interaction(map_trunc,exp), linetype=interaction(map_trunc,exp)),binwidth=800)+
        geom_vline(aes(xintercept=400), linetype="dotted")+
        ylab("Fragment count")+
        xlab("Fragment length (bp)")+
        labs(color = "Mapping status\nand experiment", linetype = "Mapping status\nand experiment")+     
        scale_linetype_manual(values=c(1,1,2,2))
#     ggsave("/hps/nobackup/goldman/ipoetzsch/boss_simulation_example/results/comp_accept_reject/2025-07-17_log_reads_accept_reject_comp_dist_by_map.png", length_dist_by_map)    

    # Plot fragment length distribution for the reads are rejected and don't map
    length_dist_rej <- joint %>%
        filter(fate==0&map_trunc==1)%>%
        ggplot()+
        geom_freqpoly(aes(length_contr, after_stat(count), group=desc(exp), colour=exp, linetype=exp), alpha=0.8,binwidth=800)+
        # scale_y_continuous(trans = 'log10')+
        geom_vline(aes(xintercept=400), linetype="dotted")+
        ylab("Fragment count")+
        xlab("Fragment length (bp)")+
        labs(color = "Experiment", linetype = "Experiment")+     
        scale_linetype_manual(name="Experiment", values=c(1,2))
    ggsave("/hps/nobackup/goldman/ipoetzsch/boss_simulation_example/results/comp_accept_reject/2025-07-17_log_reads_accept_reject_comp_dist_rej.png", 
    length_dist_rej)    

    # Plot fragment length distribution for the reads are rejected and don't map
    length_dist_rej <- joint %>%
        filter(fate==0&map_trunc==1)%>%
        ggplot()+
        geom_freqpoly(aes(length_contr, after_stat(count), group=desc(exp), colour=exp, linetype=exp), alpha=0.8,binwidth=800)+
        # scale_y_continuous(trans = 'log10')+
        geom_vline(aes(xintercept=400), linetype="dotted")+
        ylab("Fragment proportion")+
        xlab("Fragment length (bp)")+
        labs(color = "Experiment", linetype = "Experiment")+     
        scale_linetype_manual(name="Experiment", values=c(1,2))
    ggsave("/hps/nobackup/goldman/ipoetzsch/boss_simulation_example/results/comp_accept_reject/2025-07-17_log_reads_accept_reject_comp_dist_rej.png", 
    length_dist_rej)    

    # Adding time manually
    alpha <- 300
    rho <- 300
    mu <- 400
    pores <- 512
    speed <- 400
    
    # Visualise change in read number over experiment time
    time_by_n <- joint %>%
        group_by(exp)%>%
        mutate(time = length + alpha) %>%
        mutate(time_contr = length_contr + alpha)%>%
        mutate(time = if_else(length == mu, time+rho, time))%>%
        mutate(time= time/pores/speed/60, time_contr = time_contr/pores/speed/60)%>%
        mutate(n = row_number())%>%
        mutate(time = cumsum(time), time_contr = cumsum(time_contr))%>%
        # summarise(time=sum(time))%>%
        # print(width=Inf)%>%
    ggplot(data=.)+
    geom_line(aes(time, n, group=exp, colour=exp, linetype=exp), alpha=0.8)+
    geom_line(aes(time_contr, n, group=exp, colour="control", linetype="control"), alpha=0.8)+
    ylab("Fragment count")+
    xlab("Seq. time (pseudominutes)")+
    labs(color  = "Experiment", linetype = "Experiment")+
    scale_linetype_manual(name="Experiment",values=c(1,1,2))+
    # scale_linetype_manual(name="Experiment",values=c(1,2))+
    scale_color_discrete()

    ggsave(output, time_by_n)    

    # Plot sequence acquired by sequencing time
    seq <- joint %>%
        group_by(exp)%>%
        mutate(time = length + alpha) %>%
        mutate(time = if_else(length == mu, time+rho, time))%>%
        mutate(time_contr = length_contr + alpha)%>%
        mutate(time= time/pores/speed/60, time_contr = time_contr/pores/speed/60)%>%
        mutate(time = cumsum(time), length=cumsum(length), length_contr = cumsum(length_contr), time_contr=cumsum(time_contr))%>%
        # summarise(time=sum(time))%>%
        # print(width=Inf)%>%
    ggplot(data=.)+
    geom_line(aes(time_contr, length_contr, group=exp, colour="control", linetype="control"), alpha=0.8)+
    geom_line(aes(time, length, group=exp, colour=exp, linetype=exp), alpha=0.8)+
    xlab("Seq. time (pseudominutes)")+
    ylab("Total seq. read length (bp)")+
    labs(color  = "Experiment", linetype = "Experiment")+
    scale_linetype_manual(name="Experiment",values=c(1,1,2))
    # scale_linetype_manual(name="Experiment",values=c(1,2))


    # Plot sequence acquired by read number
    seq_by_n <- joint %>%
        group_by(exp)%>%
        mutate(time = length + alpha) %>%
        mutate(time = if_else(length == mu, time+rho, time))%>%
        mutate(time= time/pores/speed/60)%>%
        mutate(n = row_number())%>%
        mutate(time = cumsum(time), length=cumsum(length), length_contr=cumsum(length_contr))%>%
        # summarise(time=sum(time))%>%
        # print(width=Inf)%>%
    ggplot(data=.)+
    geom_line(aes(n, length, group=desc(exp), colour=exp, linetype=exp), alpha=0.8)+
    geom_line(aes(n, length_contr, group=desc(exp), colour="control", linetype="control"), alpha=0.8)+
    xlab("Fragment count")+
    ylab("Total seq. read length (bp)")+
    labs(color  = "Experiment", linetype = "Experiment")+
    scale_linetype_manual(name="Experiment",values=c(1,1,2))
    # scale_linetype_manual(name="Experiment",values=c(1,2))
    # ggsave(output, seq_by_n)

    # How much of the total control read length, do the two conditions read
    joint%>%
        group_by(exp)%>%
        summarise(length = sum(length), length_contr = sum(length_contr))%>%
        mutate(length_prop = length/length_contr)%>%
        print()

    
    # Output plots
    layout <- wrap_plots(list(length_dist, time_by_n, seq, seq_by_n)) +
    plot_layout(
    nrow = 2,
    ncol = 2,
    # axis_titles = "collect",
    guides = "collect",
    # axes = "collect"
    ) +
    plot_annotation(
    tag_levels = "A"
    ) &
    # scale_fill_viridis_c(limits = c(min_z, max_z)) &
    coord_cartesian(expand = FALSE)

    # layout <- (((length_dist + time_by_n) / (seq + seq_by_n)) / guide_area()) +
    # plot_annotation(tag_levels = "A") +
    # plot_layout(guides="auto") &
    # # labs(color  = "Experiment", linetype = "Experiment") &
    # theme(
    #     legend.position = "bottom",
    #    legend.title = element_blank(),
    #   strip.text.x = element_blank(),
    #   plot.background = element_rect(fill = "white", color = NA)
    # ) 

  ggsave(output, layout)
}



# input to this script is the logfile produced by the run
analyse_log()