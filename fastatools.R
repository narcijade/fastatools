#!/usr/bin/env Rscript
# fastatools implementation in R
# INPUT: FASTA file
# OUTPUT: FASTA statistics or FASTA file

###############################################################################
##### LIBRARIES ###############################################################
###############################################################################
msg.trap <- capture.output( suppressMessages(library(tidyverse)))
options(warn=-1)

###############################################################################
##### FUNCTIONS ###############################################################
###############################################################################
read_fasta <- function(filename) {
# INPUT: FASTA file
# OUTPUT: dataframe of FASTA
    raw_fasta <- read_file(filename)

    # format fasta dataframe into header, description and seq
    # description is everything after the first space in the header
    # calculate per sequence statistics
    fasta <- data.frame(all = raw_fasta) %>%
        separate_rows(all, sep=">") %>%
        filter(all != "") %>%
        separate(all, c("header", "seq"), sep = "\n", extra="merge") %>%
        mutate(seq = str_replace_all(seq, "\n", "")) %>%
        separate(header, c("seq_id", "description"), extra = "merge", sep = " ") %>%
        mutate(seq_len = nchar(seq)) %>%
        mutate(num_a = str_count(toupper(seq), "A")) %>%
        mutate(num_t = str_count(toupper(seq), "T")) %>%
        mutate(num_g = str_count(toupper(seq), "G")) %>%
        mutate(num_c = str_count(toupper(seq), "C")) %>%
        mutate(num_n = str_count(toupper(seq), "N")) %>%
        mutate(not_atgcn = str_count(toupper(seq), "[^ATGCN]")) %>%
        mutate(perc_gc = round((num_g + num_c) / (num_a + num_t + num_g + num_c), 2)) %>%
        mutate(perc_n = round(num_n / seq_len, 2))
}


# write out fasta file after filtering
write_fasta <- function (df) {
# INPUT: FASTA dataframe
# OUTPUT: FASTA to STDOUT
    out <- df %>%
        select(seq_id, description, seq) %>%
        mutate(header = ifelse(is.na(description), paste0(">", seq_id) ,paste0(">", seq_id, " ", description))) %>%
        select(header, seq) %>%
        mutate(seq = str_replace_all(seq, "(.{60})", "\\1\n")) %>%
        mutate(line_out = paste0(header, "\n", seq, "\n")) %>%
        select(line_out)

    return(cat(out$line_out, sep=''))
}

stats <- function(filename) {
# Calculate summary statistics
    fasta <- read_fasta(filename)
    total_length <- sum(fasta$seq_len)
    ordered_lengths <- rev(sort(fasta$seq_len))
    n50_pos <- cumsum(ordered_lengths) <= total_length / 2
    n50 <- ordered_lengths[n50_pos]

    fasta_stats <- fasta %>%
        summarise(perc_gc = round((sum(num_g) + sum(num_c)) / (sum(num_a) + sum(num_t) + sum(num_g) + sum(num_c)), 2),
                  perc_n = round(sum(num_n) / total_length, 2),
                  longest = max(seq_len),
                  shortest = min(seq_len),
                  average = mean(seq_len)
                  )

    # Write stats
    cat( "Total Length: ", total_length, "\n",
         " Mean Length: ", fasta_stats$average,"\n",
         " Longest Seq: ", fasta_stats$longest,"\n",
         "Shortest Seq: ", fasta_stats$shortest,"\n",
         "         %GC: ", fasta_stats$perc_gc * 100,"\n",
         "          %N: ", fasta_stats$perc_n * 100, "\n",
         "         N50: ", n50, "\n", sep = ''
    )
}

gc_hist <- function(filename) {
# GC histogram
    msg.trap <- capture.output( suppressMessages(library(cowplot)))
    msg.trap <- capture.output( suppressMessages(library(gridExtra)))
    msg.trap <- capture.output( suppressMessages(library(R.devices)))

    fasta <- read_fasta(filename)
    theme_set(theme_cowplot())
    p <- ggplot(fasta, aes(x = perc_gc)) +
        geom_histogram(binwidth=.01) +
        xlim(-.01,1.01) +
        xlab("%GC") +
        ylab("# of Seqs") +
        ggtitle("Histogram of GC content")
    # write histogram out
    png_name <- str_replace(filename, "(.*)\\..*", "\\1_gc_hist.png")
    suppressGraphics(ggsave(png_name, plot = p, width = 5, height = 5, units = "in"))
}

len_hist <- function(filename) {
# Length histogram
    msg.trap <- capture.output( suppressMessages(library(cowplot)))
    msg.trap <- capture.output( suppressMessages(library(gridExtra)))
    msg.trap <- capture.output( suppressMessages(library(R.devices)))

    fasta <- read_fasta(filename)
    theme_set(theme_cowplot())
    p <- ggplot(fasta, aes(x = seq_len)) +
        geom_histogram(binwidth=100) +
        xlab("Seq Length") +
        ylab("# of Seqs") +
        ggtitle("Histogram of Sequence Length")
    # write histogram out
    png_name <- str_replace(filename, "(.*)\\..*", "\\1_len_hist.png")
    suppressGraphics(ggsave(png_name, plot = p, width = 5, height = 5, units = "in"))
}

len_filter <- function(filename, len_low, len_high) {
# Filter FASTA by length
    fasta <- read_fasta(filename)

    filtered_by_length <- fasta %>%
        filter(seq_len <= len_high) %>%
        filter(seq_len >= len_low)

    write_fasta(filtered_by_length)
}

gc_filter <- function(filename, gc_low, gc_high) {
# Filter FASTA by gc
    fasta <- read_fasta(filename)
    filtered_by_gc <- fasta %>%
        filter(perc_gc <= gc_high) %>%
        filter(perc_gc >= gc_low)

    write_fasta(filtered_by_gc)
}

###############################################################################
##### MAIN ####################################################################
###############################################################################

args = commandArgs(trailingOnly=TRUE)

if (args[1] == "stats" & length(args) == 2) {
    stats(args[2])
} else if (args[1] == "gc_hist" & length(args) == 2) {
    gc_hist(args[2])
} else if (args[1] == "len_hist" & length(args) == 2) {
    len_hist(args[2])
} else if (args[1] == "gc_filter" &
           length(args) == 4 &
           is.numeric(args[3]) &
           is.numeric(args[4])
          ) {
    gc_filter(args[2], args[3], args[4])
} else if (args[1] == "len_filter" &
           length(args) == 4 &
           is.numeric(args[3]) &
           is.numeric(args[4])
          ) {
    len_filter(args[2], args[3], args[4])
} else {
    stop("\nUsage:\nstats FASTA\ngc_hist FASTA\nlen_hist FASTA\ngc_filter FASTA min_gc max_gc\nlen_filter min_len max_len\n", call.=FALSE)
}
