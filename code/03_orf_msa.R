library(Biostrings)
library(DECIPHER)
library(dplyr)
library(data.table)
library(ggplot2)
library(BuenColors)

# note: the csv has its own `genus` column, so keep this map named differently
gmap <- c(MZ286102 = "alpha", MZ286188 = "alpha",
          MZ286258 = "beta",  OZ258728 = "beta",
          MZ286197 = "gamma", OP549859 = "gamma")

cols <- c(alpha = "#984EA3", beta = "#7A0403FF", gamma = "#4777EFFF",
          mismatch = "grey86", gap = "white")

ss <- fread("../data/selected_strains.csv") %>%
  mutate(virus = gsub("\\.1$", "", target_id), gen = unname(gmap[virus])) %>%
  arrange(gen, virus)

align_1 <- function(seqs, nms){
  aa <- AAStringSet(seqs); names(aa) <- nms
  m <- do.call(rbind, strsplit(as.character(AlignSeqs(aa, verbose = FALSE)), ""))
  rownames(m) <- nms

  cons <- apply(m, 2, function(v){
    v <- v[v != "-"]
    if (!length(v)) return("-")
    names(sort(table(v), decreasing = TRUE))[1]
  })

  data.frame(virus = rep(rownames(m), ncol(m)),
             pos   = rep(seq_len(ncol(m)), each = nrow(m)),
             aa    = as.vector(m)) %>%
    mutate(fill = case_when(aa == "-" ~ "gap",
                            aa == cons[pos] ~ unname(gmap[virus]),
                            TRUE ~ "mismatch"))
}

msa_df <- rbind(align_1(ss$orf1_seq, ss$virus) %>% mutate(orf = "ORF1"),
                align_1(ss$orf2_seq, ss$virus) %>% mutate(orf = "ORF2")) %>%
  mutate(virus = factor(virus, levels = rev(ss$virus)),
         fill  = factor(fill, levels = c("alpha", "beta", "gamma", "mismatch", "gap")))

msa_1 <- function(which_orf, legend = FALSE){
  d <- msa_df %>% filter(orf == which_orf)
  ycol <- cols[gmap[levels(d$virus)]]

  p <- ggplot(d, aes(x = pos, y = virus, fill = fill)) +
    geom_tile() +
    scale_fill_manual(values = cols, drop = FALSE,
                      labels = c("alpha match", "beta match", "gamma match", "mismatch", "gap")) +
    scale_x_continuous(expand = c(0, 0)) +  scale_y_discrete(expand = c(0, 0)) +
    labs(x = "", y = "", title = paste0(which_orf, " (", max(d$pos), " aa)")) +
    pretty_plot(fontsize = 6) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6, face = "plain", margin = margin(b = 1)),
          plot.margin = margin(2, 4, 0, 1),
          axis.text.y = element_text(size = 5, color = ycol),
          axis.text.x = element_text(size = 5),
          axis.ticks.y = element_blank())

  if (!legend) return(p)

  p + guides(fill = guide_legend(nrow = 1, title = NULL, keywidth = unit(6, "pt"),
                                 keyheight = unit(5, "pt"))) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 5),
          legend.key = element_rect(fill = NA, colour = "grey40", linewidth = 0.2),
          legend.margin = margin(0, 0, 0, 0))
}

p <- cowplot::plot_grid(msa_1("ORF1"), msa_1("ORF2"), ncol = 1, align = "v")
p
cowplot::ggsave2(p, file = "../plots/orf_msa.pdf", width = 3.7, height = 1.8)

leg <- cowplot::get_plot_component(msa_1("ORF2", legend = TRUE), "guide-box-bottom")
p_leg <- cowplot::plot_grid(msa_1("ORF1"), msa_1("ORF2"), leg, ncol = 1,
                            rel_heights = c(1, 1, 0.2))
cowplot::ggsave2(p_leg, file = "../plots/orf_msa_legend.pdf", width = 4.3, height = 2.1)

msa_df %>% group_by(orf, virus) %>%
  summarize(pct_consensus = round(100 * mean(!fill %in% c("gap", "mismatch")), 1),
            pct_gap = round(100 * mean(fill == "gap"), 1), .groups = "drop") %>%
  data.frame()

library(Biostrings)
library(pwalign)
library(dplyr)
library(data.table)

gmap <- c(MZ286102 = "alpha", MZ286188 = "alpha",
          MZ286258 = "beta",  OZ258728 = "beta",
          MZ286197 = "gamma", OP549859 = "gamma")

ss <- fread("../data/selected_strains.csv") %>%
  mutate(virus = gsub("\\.1$", "", target_id), gen = unname(gmap[virus])) %>%
  arrange(gen, virus)

data(BLOSUM62)

pw <- function(seqs, nms){
  s   <- setNames(AAStringSet(seqs), nms)
  len <- setNames(nchar(seqs), nms)                 # plain integer vector
  g   <- expand.grid(a = nms, b = nms, stringsAsFactors = FALSE) %>% filter(a < b)
  
  lapply(seq_len(nrow(g)), function(i){
    al <- pairwiseAlignment(s[[g$a[i]]], s[[g$b[i]]], substitutionMatrix = BLOSUM62,
                            gapOpening = 10, gapExtension = 0.5, type = "global")
    pc <- strsplit(as.character(alignedPattern(al)), "")[[1]]
    qc <- strsplit(as.character(alignedSubject(al)), "")[[1]]
    ung <- pc != "-" & qc != "-"
    
    data.table(a = g$a[i], b = g$b[i],
               pid_aln   = 100 * sum(pc == qc & ung) / length(pc),
               pid_short = 100 * sum(pc == qc & ung) / min(len[g$a[i]], len[g$b[i]]),
               sim_pct   = 100 * mean(mapply(function(x, y) BLOSUM62[x, y] > 0, pc[ung], qc[ung])),
               aln_len   = length(pc))
  }) %>% rbindlist()
}

lapply(c("orf1_seq", "orf2_seq"), function(o){
  pw(ss[[o]], ss$virus) %>%
    mutate(orf = toupper(sub("_seq", "", o)),
           pair = ifelse(gmap[a] == gmap[b], paste0("within-", gmap[a]), "between"))
}) %>% rbindlist() -> pid_df

pid_df %>% filter(pair != "between") %>% arrange(orf, pair) %>% data.frame()
pid_df %>% group_by(orf, pair) %>% summarize(median_pid = median(pid_short), .groups = "drop")