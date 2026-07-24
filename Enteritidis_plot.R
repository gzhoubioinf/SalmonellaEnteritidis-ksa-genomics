library(tidyverse)
library(ape)
library(ggtree)
library(ggnewscale)
library(rhierbaps)
library(phytools)
library(phangorn)
library(cowplot)
library(patchwork)

# Figure1A
metadata_hierbaps <- read_tsv("./Figure1/metadata_baps.tsv", show_col_types = FALSE)
annot <- metadata_hierbaps %>%
  transmute(
    Isolate,
    Origin = factor(str_to_title(origin), levels = c("External", "Internal")),
    Year   = factor(as.numeric(format(as.Date(collection_date, "%d/%m/%Y"), "%Y")),
                    levels = c(2009, 2014:2024)),
    Group  = factor(baps_level2, levels = as.character(1:12)),
    Source = factor(ifelse(is.na(Source) | Source == "NA", "NA", Source),
                    levels = c("Human", "Food", "Environment", "NA")),
    Region = factor(region)
  ) %>%
  column_to_rownames("Isolate")

col_origin <- c(External = "#4299E1", Internal = "#a6cee3")
col_source <- c("Human"="blue", "Food"="red", "Environment"="#FFD700", "NA"="grey70")
col_group  <- c("1"="#5d2a4f","2"="#8c510a","3"="#b35fb0","4"="#5aa54a",
                "5"="#4f51a0","6"="#f1a14a","7"="#fff200","8"="#e8467f",
                "9"="#349a93","10"="#b6b3e3","11"="#a13a36","12"="#74ef3f")
col_region <- c("Africa"="#1f78b4","East Asia"="#33a02c","Europe"="#6a3d9a",
                "Middle East"="#1b9e77","North America"="#ff7f00","Oceania"="#a6761d")
col_year   <- c("2009"="#5e4fa2","2014"="#3288bd","2015"="#66c2a5","2016"="#abdda4",
                "2017"="#e6f598","2018"="#ffffbf","2019"="#fee08b","2020"="#fdae61",
                "2021"="#f46d43","2022"="#d53e4f","2023"="#9e0142","2024"="#67001f")

dna    <- read.dna("Figure1/Enteritidis_snp_sites.aln", format = "fasta")
njtree <- midpoint.root(nj(dist.dna(dna, model = "N", pairwise.deletion = TRUE)))
p <- ggtree(njtree, size = 0.4) +
  geom_treescale(width = 4, x = 0, y = Ntip(njtree) * 0.75, offset = 2, fontsize = 4)

tw <- max(p$data$x, na.rm = TRUE)
sw <- tw * 0.14; gap <- tw * 0.006; step <- sw + gap
off <- function(i) tw * 0.005 + (i - 1) * step
add_strip <- function(g, col, i)
  gheatmap(g, annot[col], offset = off(i), width = sw / tw, color = NA, colnames = FALSE)

p1 <- add_strip(p,                "Origin", 1) + scale_fill_manual(values = col_origin,
                                                                   na.value = "grey80")
p2 <- add_strip(p1 + new_scale_fill(), "Year",   2) + scale_fill_manual(values = col_year,
                                                                        na.value = "grey80")
p3 <- add_strip(p2 + new_scale_fill(), "Region", 3) + scale_fill_manual(values =
                                                                          col_region, na.value = "grey80")
p4 <- add_strip(p3 + new_scale_fill(), "Source", 4) + scale_fill_manual(values =
                                                                          col_source, na.value = "grey80")
p5 <- add_strip(p4 + new_scale_fill(), "Group",  5) + scale_fill_manual(values = col_group,
                                                                        na.value = "grey80")

strip_names <- c("Origin", "Year", "Region", "Source", "Group")
strip_x <- sapply(seq_along(strip_names), function(i) tw + off(i) + sw / 2)
p_labeled <- p5 +
  annotate("text", x = strip_x + tw * 0.06, y = -2, label = strip_names,
           angle = 45, hjust = 1, size = 3, colour = "grey20",
           family = "sans", fontface = "bold")

legend_theme <- theme(
  legend.justification = c(0, 1),
  legend.key.size = unit(0.35, "cm"),
  legend.title    = element_text(size = 8, face = "bold"),
  legend.text     = element_text(size = 7),
  legend.margin   = margin(0, 0, 0, 0)
)

make_legend <- function(values, breaks, title, na = FALSE) {
  if (na) { values <- c(values, "NA" = "grey80"); breaks <- c(breaks, "NA") }
  df <- data.frame(value = factor(breaks, levels = breaks))
  get_legend(
    ggplot(df, aes(1, value, fill = value)) + geom_tile() +
      scale_fill_manual(values = values, breaks = breaks, name = title) +
      theme_void() + legend_theme
  )
}

legend_origin <- make_legend(col_origin, c("External", "Internal"), "Origin")
legend_year   <- make_legend(col_year,   names(col_year),  "Year",   na = TRUE)
legend_region <- make_legend(col_region, names(col_region),"Region")
legend_source <- make_legend(col_source, c("Human","Food","Environment","NA"), "Source")
legend_group  <- make_legend(col_group,  as.character(1:12), "Group")

col_left  <- plot_grid(legend_origin, legend_year, legend_region,
                       ncol = 1, align = "v", axis = "l", rel_heights = c(3, 8, 4))
col_right <- plot_grid(legend_source, legend_group, NULL,
                       ncol = 1, align = "v", axis = "l", rel_heights = c(5, 16, 4))
legend_block <- plot_grid(col_left, col_right, nrow = 1, rel_widths = c(0.85, 1))
legend_col   <- plot_grid(legend_block, NULL, ncol = 1, rel_heights = c(1, 0.15))

p_main <- p_labeled +
  theme(legend.position = "none", plot.margin = margin(5, 0, 30, 5)) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06)))

fig <- plot_grid(p_main, legend_col, nrow = 1, rel_widths = c(3, 1.4),
                  label_size = 16, label_fontface = "bold",
                 label_x = 0, label_y = 1)

library(adegenet)

snp <- suppressWarnings(
  fasta2genlight("Figure1/Enteritidis_snp_sites.aln", snpOnly = TRUE,
                 chunkSize = 100, quiet = TRUE)
)
stopifnot(nInd(snp) == 397, nLoc(snp) == 1218)
pca <- glPca(snp, nf = 10)

ve   <- 100 * pca$eig / sum(pca$eig)
xlab <- paste0("PC1 (", round(ve[1], 1), "%)")
ylab <- paste0("PC2 (", round(ve[2], 1), "%)")

pca_df <- as_tibble(pca$scores[, 1:2], rownames = "Isolate") %>%
  left_join(metadata_hierbaps, by = "Isolate") %>%
  mutate(Source = factor(Source, levels = c("Human", "Food", "Environment")),
         Group  = factor(baps_level2, levels = as.character(1:12)))

col_source_bc <- c(Human = "blue", Food = "red", Environment = "#FFD700")

pB <- ggplot(arrange(pca_df, Source), aes(PC1, PC2, colour = Source)) +
  geom_point(size = 2, alpha = 0.7,
             position = position_jitter(width = 0.04, height = 0.04, seed = 1)) +
  scale_colour_manual(values = col_source_bc, na.translate = FALSE) +
  labs(x = xlab, y = ylab, title = "B") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.title.position = "plot", legend.position = "none")

zoom_bc <- pca_df %>% filter(Group %in% as.character(1:8))
pad     <- 0.05
xlim_c  <- range(zoom_bc$PC1) + c(-1, 1) * pad * diff(range(zoom_bc$PC1))
ylim_c  <- range(zoom_bc$PC2) + c(-1, 1) * pad * diff(range(zoom_bc$PC2))

pC_base <- ggplot(pca_df, aes(PC1, PC2, colour = Group)) +
  geom_point(size = 2, alpha = 0.7,
             position = position_jitter(width = 0.04, height = 0.04, seed = 1)) +
  scale_colour_manual(values = col_group, drop = FALSE) +
  theme_classic() +
  theme(legend.position = "none")

pC_overview <- ggplot(pca_df, aes(PC1, PC2, colour = Group)) +
  geom_point(size = 0.7, alpha = 0.8,
             position = position_jitter(width = 0.04, height = 0.04, seed = 1)) +
  scale_colour_manual(values = col_group, drop = FALSE) +
  theme_classic() +
  theme(legend.position = "none") +
  annotate("rect", xmin = xlim_c[1], xmax = xlim_c[2],
           ymin = ylim_c[1], ymax = ylim_c[2],
           fill = NA, colour = "grey30", linewidth = 0.3) +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.4),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(1, 1, 1, 1))

pC <- pC_base +
  coord_cartesian(xlim = xlim_c, ylim = ylim_c) +
  labs(x = xlab, y = ylab, title = "C") +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.title.position = "plot") +
  inset_element(pC_overview, left = 0.00, bottom = 0.02, right = 0.36, top = 0.39,
                align_to = "panel")

panel_BC <- plot_grid(pB, pC, ncol = 1, rel_heights = c(1, 1))
Figure1 <- plot_grid(fig, panel_BC, nrow = 1, rel_widths = c(2, 1))

Figure1 <- ggdraw(Figure1) +
  draw_plot_label("A", x = 0.005, y = 0.99, size = 16, fontface = "bold",
                  hjust = 0, vjust = 1)

ggsave("./Figure1/Figure1.svg", Figure1, width = 8.67, height = 4.6, device = svglite::svglite)

#Figure2
library(ape); library(dplyr); library(tidyr); library(readr)
library(ggtree); library(phytools); library(ggplot2); library(aplot)
library(RColorBrewer)

star <- function(p) ifelse(is.na(p), "",
                      ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                      ifelse(p < 0.05, "*", ""))))

vs_rest <- function(df, var, group = "BAPS") {
  g <- as.character(df[[group]]); x <- df[[var]]
  res <- lapply(sort(unique(g)), function(k) {
    a <- x[g == k]; b <- x[g != k]
    p <- if (length(na.omit(a)) < 1 || length(unique(na.omit(c(a, b)))) < 2) NA_real_
         else wilcox.test(a, b, alternative = "greater", exact = FALSE)$p.value
    data.frame(BAPS = k, median_group = median(a, na.rm = TRUE),
               median_rest = median(b, na.rm = TRUE), p = p)
  })
  out <- do.call(rbind, res)
  out$q   <- p.adjust(out$p, method = "BH")
  out$lab <- star(out$q)
  out$xpos <- max(x, na.rm = TRUE) * 1.02
  out
}

  dna     <- read.dna("Figure1/Enteritidis_snp_sites.aln", format = "fasta")
  distdna <- dist.dna(dna, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
  rownames(distdna) <- as.character(rownames(distdna))
  colnames(distdna) <- rownames(distdna)

  meta <- read_tsv("Figure1/metadata_baps.tsv", show_col_types = FALSE) %>%
    mutate(BAPS = as.character(baps_level2))

  region_lv  <- c("Middle East", "Europe", "North America", "East Asia",
                  "Oceania", "Africa")
  country_lv <- c("Saudi Arabia", "Kuwait",
                  "United Kingdom",
                  "United States",
                  "China",
                  "Other")
  col_region  <- c("Africa"="#1f78b4","East Asia"="#33a02c","Europe"="#6a3d9a",
                   "Middle East"="#1b9e77","North America"="#ff7f00","Oceania"="#a6761d")
  col_region  <- col_region[region_lv]
  col_source  <- c("Human"="blue","Food"="red","Environment"="#FFD700","NA"="grey70")
  col_country <- c(
    "China"          = "#e78ac3",
    "Kuwait"         = "#a6d854",
    "Saudi Arabia"   = "#8da0cb",
    "United Kingdom" = "#fc8d62",
    "United States"  = "#66c2a5",
    "Other"          = "#b3b3b3"
  )
  col_country <- col_country[country_lv]

  meta <- meta %>%
    mutate(Country = factor(ifelse(location %in% names(col_country), location, "Other"),
                            levels = names(col_country)))

  meta <- meta %>%
    mutate(Source = ifelse(is.na(Source) | Source == "NA", "NA", Source),
           Source = factor(Source, levels = c("Human", "Food", "Environment", "NA")),
           region = factor(region, levels = region_lv))

  top_title_bold <- function(t) list(ggtitle(t),
                                     theme(plot.title = element_text(hjust = 0.5, size = 7, face =
                                                                       "bold")))

  big_key <- guides(fill = guide_legend(override.aes = list(size = 4)))

  bar_theme <- theme_bw() +
    theme(axis.text.x = element_text(size = 8, colour = NA),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(), axis.title =
            element_blank(),
          legend.position = "right",
          legend.title    = element_text(size = 7, face = "bold"),
          legend.text     = element_text(size = 7),
          legend.key.size = unit(0.22, "cm"),
          legend.spacing.y = unit(0.2, "cm"),
          plot.margin = margin(8, 2, 8, 0))

  strip_theme <- theme_bw() +
    theme(axis.text.x  = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.text.y  = element_blank(), axis.title.y = element_blank(),
          plot.margin  = margin(8, 10, 8, 10))

  clusters <- as.character(sort(unique(as.numeric(meta$BAPS))))
  avg_d <- function(a, b) {
    s1 <- meta$Isolate[meta$BAPS == a]; s2 <- meta$Isolate[meta$BAPS == b]
    mean(distdna[as.character(s1), as.character(s2)], na.rm = TRUE)
  }
  M <- matrix(0, length(clusters), length(clusters), dimnames = list(clusters, clusters))
  for (i in clusters) for (j in clusters) if (i != j) M[i, j] <- avg_d(i, j)
  M <- apply(M, c(1, 2), round)

  njtree_baps <- midpoint.root(bionj(M))
  p1 <- ggtree(njtree_baps) + geom_treescale() +
    geom_tiplab(aes(label = label), align = TRUE, size = 2.8, offset = 0.2, linesize = .2) +
    xlim(0, max(node.depth.edgelength(njtree_baps)) * 1.1) +
    theme(plot.margin = margin(8, 0, 8, 8))

  pRegion  <- ggplot(meta, aes(1, BAPS, fill = region)) +
    geom_col(position = position_fill(reverse = TRUE), width = 1) +
    scale_fill_manual(name = "Region", values = col_region) + bar_theme + big_key +
    top_title_bold("Region")

  pCountry <- ggplot(meta, aes(1, BAPS, fill = Country)) +
    geom_col(position = position_fill(reverse = TRUE), width = 1) +
    scale_fill_manual(name = "Country", values = col_country) + bar_theme + big_key +
    top_title_bold("Country")

  pSource  <- ggplot(meta, aes(1, BAPS, fill = Source)) +
    geom_col(position = position_fill(reverse = TRUE), width = 1) +
    scale_fill_manual(name = "Source", values = col_source) + bar_theme + big_key +
    top_title_bold("Source")

  cnt <- meta %>% count(BAPS)
  p2 <- ggplot(cnt, aes(n, BAPS)) + geom_col(fill = "blue") +
    strip_theme + top_title_bold("Sample count")

  age <- meta %>%
    separate(collection_date, c("d", "m", "y"), sep = "/", remove = TRUE) %>%
    mutate(age = 2025 - as.numeric(y))
  age_st <- vs_rest(age, "age")
  p4 <- ggplot(age, aes(age, BAPS)) + geom_boxplot(fill = "#00CC66") +
    geom_text(data = age_st, aes(x = xpos, y = BAPS, label = lab),
              hjust = 0, size = 2.8, inherit.aes = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    strip_theme + top_title_bold("Age (years)")

  dist_df <- as.data.frame(as.table(distdna)) %>%
    rename(s1 = Var1, s2 = Var2, snp = Freq) %>%
    mutate(s1 = as.character(s1), s2 = as.character(s2)) %>%
    filter(s1 < s2)
  tr <- meta %>% select(Isolate, BAPS)
  within <- dist_df %>%
    left_join(tr, by = c("s1" = "Isolate")) %>% rename(B1 = BAPS) %>%
    left_join(tr, by = c("s2" = "Isolate")) %>% rename(B2 = BAPS) %>%
    filter(B1 == B2) %>% rename(BAPS = B1)
  snp_st <- vs_rest(within, "snp")
  p6 <- ggplot(within, aes(snp, BAPS)) +
    geom_boxplot(outlier.shape = NA, width = .6, fill = "grey80",
                 linewidth = 0.25, fatten = 1) +
    geom_jitter(height = .15, width = 0, alpha = .05, size = .3) +
    geom_text(data = snp_st, aes(x = xpos, y = BAPS, label = lab),
              hjust = 0, size = 2.8, inherit.aes = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    strip_theme + top_title_bold("Pairwise SNP Distance")

  amr <- read_tsv("Figure1/amrfinder_all.tsv", show_col_types = FALSE) %>%
    filter(`Gene symbol` != "Gene symbol") %>%
    rename(Gene_symbol = `Gene symbol`, Element_type = `Element type`)
  amr_count <- amr %>% filter(Element_type == "AMR") %>%
    distinct(Sample_name, Gene_symbol) %>%
    count(Sample_name, name = "AMR_sum")
  baps_amr <- meta %>%
    left_join(amr_count, by = c("Isolate" = "Sample_name")) %>%
    mutate(AMR_sum = replace_na(AMR_sum, 0))
  amr_st <- vs_rest(baps_amr, "AMR_sum")
  p7 <- ggplot(baps_amr, aes(AMR_sum, BAPS)) + geom_boxplot(fill = "#107E7D") +
    geom_text(data = amr_st, aes(x = xpos, y = BAPS, label = lab),
              hjust = 0, size = 2.8, inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq(0, 12, 2),
                       limits = c(NA, max(baps_amr$AMR_sum) * 1.12)) +
    strip_theme + top_title_bold("AMR count")

  fig2 <- p7 %>%
    insert_left(p6, width = 1.5) %>%
    insert_left(p4) %>%
    insert_left(p2) %>%
    insert_left(pSource,  width = 0.3) %>%
    insert_left(pCountry, width = 0.3) %>%
    insert_left(pRegion,  width = 0.3) %>%
    insert_left(p1, width = 1.2)

  fig2
  ggsave("Figure1/Figure2.svg", fig2, width = 10, height = 4.1, device = svglite::svglite)

  #Figure3

  library(readxl)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggnewscale)
  options(scipen = 999)
  #Figure3A
  beast <- read_excel("./Figure1/BAPS_Age_clock_results.xlsx") %>%
    mutate(BAPS = factor(BAPS, levels = sort(unique(as.numeric(BAPS)))))

  big_theme <- theme_bw() +
    theme(text       = element_text(size = 10, face = "bold",colour = "black"),
          axis.text  = element_text(size = 10,face = "plain",colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold", colour = "black"),
          plot.margin = margin(10, 25, 10, 20))

  pAge <- ggplot(beast, aes(BAPS, age)) +
    geom_bar(stat = "identity", fill = "#77c3e5") +
    geom_errorbar(aes(ymin = age_lower, ymax = age_upper),
                  width = 0.4, colour = "black", linewidth = 0.5, alpha = 0.9) +
    coord_cartesian(ylim = c(2005, 2024)) +
    scale_y_continuous(breaks = seq(1990, 2024, 5), expand = c(0, 0)) +
    labs(x = "BAPS", y = "MRCA Age (year)") +
    big_theme

  pRate <- ggplot(beast, aes(BAPS, clockrate_normalized)) +
    geom_bar(stat = "identity", fill = "#77c3e5") +
    geom_errorbar(aes(ymin = clockrate_lower_normalized, ymax = clockrate_upper_normalized),
                  width = 0.4, colour = "black", linewidth = 0.5, alpha = 0.9) +
    labs(x = "BAPS", y = "Substitution Rate (subs/site/year)") +
    big_theme

  fig3a <- (pAge + pRate) +
    plot_annotation(title = "A",
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust =
                                                              0))) &
    theme(plot.margin = margin(3, 3, 3, 3))
  fig3a
  ggsave("Figure1/Figure_3A.svg", fig3a, width = 10.67, height = 6.6, device = svglite::svglite)

  #Figure3B

  library(skygrowth)
  library(treeio)
  library(ggplot2)
  library(patchwork)
  options(scipen = 999)
  sky_theme <- theme_bw() +
    theme(axis.text    = element_text(size = 9, hjust = 1, colour = "black"),
          axis.title   = element_text(colour = "black", size = 10, face = "bold"),
          strip.text.x = element_text(colour = "black", size = 10, face = "bold"),
          plot.title   = element_text(size = 10, face = "bold", hjust = 0.5))

  no_sci <- function(x) format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
  sky_plot <- function(baps) {
    tree <- read.nexus(sprintf("Figure1/baps%s_country_beast1.tree", baps))
    fit  <- skygrowth.map(tree, res = 24, tau0 = 0.05)
    plot(fit) + sky_theme +
      xlab("Time after MRCA") +
      ggtitle(paste("BAPS", baps)) +
      scale_y_log10(labels = no_sci)
  }

  baps_groups <- c(1, 2, 3, 6, 7)
  plots <- lapply(baps_groups, sky_plot)

  combined_plot <- (plots[[1]] | plots[[2]]) /
    (plots[[3]] | plots[[4]] | plots[[5]])

  fig3b <- combined_plot +
    plot_annotation(title = "B",
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust =
                                                              0)))

  ggsave("Figure1/Figure_3B.svg", fig3b, width = 10.67, height = 6.6, device =
           svglite::svglite)

  #Figure3c
  library(ggtree); library(ape); library(treeio)
  library(patchwork); library(ggplot2); library(dplyr)
  library(readr);library(cowplot)

  all_meta <- read_tsv("Figure1/metadata_baps.tsv", show_col_types = FALSE)
  all_meta <- all_meta %>%
    mutate(sample_id = Isolate,
           Source = ifelse(is.na(Source) | Source == "NA", "NA", Source))

  country_cols <- c(
    "China"          = "#e78ac3",
    "Kuwait"         = "#a6d854",
    "Saudi_Arabia"   = "#8da0cb",
    "United_Kingdom" = "#fc8d62",
    "USA"            = "#66c2a5",
    "Netherlands"    = "#e7298a",
    "Australia"      = "#f02927",
    "Canada"         = "#1b9e77",
    "South_Africa"   = "#a9d4e7",
    "Lebanon"        = "#ffff33",
    "Iraq"           = "#b15928"
  )

  source_cols <- c("Human"="#0000ff", "Food"="#ff0000", "Environment"="#FFD700",
                   "NA"="grey70")

  country_breaks <- names(country_cols)
  country_labels <- dplyr::recode(gsub("_", " ", country_breaks), "USA" = "United States")

  make_panel <- function(baps, strip_w = 0.06) {
    tree <- read.beast(sprintf("Figure1/baps%s_country_beast1.tree", baps))
    if ("Country" %in% colnames(tree@data)) tree@data <- dplyr::rename(tree@data, country =
                                                                         Country)

    tl   <- tree@phylo$tip.label
    mrsd <- max(as.Date(sapply(strsplit(tl, "\\|"), `[`, 2)), na.rm = TRUE)
    tw   <- max(node.depth.edgelength(tree@phylo))
    mrsd_yr <- as.numeric(format(mrsd, "%Y")) + (as.numeric(format(mrsd, "%j")) - 1) / 365
    root_yr <- mrsd_yr - tw
    step    <- max(1, round(tw / 3))
    end_yr  <- floor(mrsd_yr)
    brks    <- seq(ceiling(root_yr), end_yr, by = step)
    brks    <- c(brks[end_yr - brks >= step / 2], end_yr)

    src <- data.frame(label = tl, sample_id = sub("\\|.*", "", tl), stringsAsFactors = FALSE)%>%
      left_join(all_meta %>% select(sample_id, Source), by = "sample_id")

    p <- ggtree(tree, aes(color = country), mrsd = mrsd, right = TRUE, size = 0.6) +
      theme_tree2()

    tip_df <- p$data %>%
      dplyr::filter(isTip) %>%
      dplyr::mutate(xend = max(x, na.rm = TRUE) + 0.3)

    p <- p +
      geom_segment(
        data = tip_df,
        aes(x = x, xend = xend, y = y, yend = y),
        inherit.aes = FALSE,
        linetype = "dashed",
        linewidth = 0.2,
        color = "grey60"
      ) +
      geom_range(range = "height_0.95_HPD", color = "grey30", alpha = 0.4, size = 0.8) +
      scale_color_manual(values = country_cols, na.value = "grey40",
                         breaks = country_breaks, labels = country_labels) +
      theme(axis.text.x = element_text(size = 10, colour = "black"))

    hd <- src[, c("label", "Source")]; rownames(hd) <- hd$label; hd$label <- NULL

    gheatmap(p, hd, offset = 0, width = strip_w, colnames = FALSE) +
      scale_x_continuous(breaks = brks) +
      scale_fill_manual(breaks = c("Human", "Food", "Environment"), values = source_cols) +
      ggtitle(paste("BAPS", baps)) +
      theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
            legend.position = "none")
  }

  p1 <- make_panel(1, strip_w = 0.04)
  p2 <- make_panel(2, strip_w = 0.04)
  p3 <- make_panel(3, strip_w = 0.06)
  p6 <- make_panel(6, strip_w = 0.06)
  p7 <- make_panel(7, strip_w = 0.06)

  leg_theme <- theme(legend.title = element_text(size = 11, face = "bold"),
                     legend.text  = element_text(size = 10),
                     legend.key.size = unit(0.4, "cm"),
                     legend.spacing.y = unit(1.5, "cm"))

  combo <- ggplot() +
    geom_tile(data = data.frame(v = factor(country_breaks, levels = country_breaks)),
              aes(1, v, fill = v)) +
    scale_fill_manual(values = country_cols, breaks = country_breaks,
                      labels = country_labels, name = "Country") +
    new_scale_fill() +
    geom_tile(data = data.frame(v = factor(c("Human","Food","Environment"),
                                           levels = c("Human","Food","Environment"))),
              aes(2, v, fill = v)) +
    scale_fill_manual(values = source_cols, name = "Source") +
    theme_void() + leg_theme

  legends <- get_legend(combo)

  trees <- (p1 + p2) / (p3 + p6 + p7)
  fig3c <- plot_grid(trees, legends, ncol = 2, rel_widths = c(1, 0.16))

  fig3c <- ggdraw(fig3c) +
    draw_plot_label("C", x = 0.005, y = 0.995, size = 20, fontface = "bold", hjust = 0, vjust
                    = 1)

  ggsave("Figure1/Figure_3C.svg", fig3c, width = 12, height = 10, device = svglite::svglite)

  left <- plot_grid(fig3a, fig3b, ncol = 1, rel_heights = c(0.85, 1.15))

  Figure3 <- plot_grid(left, fig3c, ncol = 2, rel_widths = c(1, 1.2))

  ggsave("Figure1/Figure3.svg", Figure3, width = 15, height = 8.0, device = svglite::svglite)

  #Figure4

  library(adegenet); library(igraph); library(ape)
  library(dplyr); library(tidyr); library(readr)
  library(rnaturalearth); library(sf)

  all_meta <- read_tsv("Figure1/metadata_baps.tsv", show_col_types = FALSE)
  all_meta$label           <- all_meta$Isolate
  all_meta$collection_date <- as.Date(all_meta$collection_date, format = "%d/%m/%Y")

  location_cols <- c(
    "China"          = "#e78ac3",
    "Kuwait"         = "#a6d854",
    "Saudi Arabia"   = "#8da0cb",
    "United Kingdom" = "#fc8d62",
    "United States"  = "#66c2a5",
    "Netherlands"    = "#e7298a",
    "Australia"      = "#f02927",
    "Canada"         = "#1b9e77",
    "South Africa"   = "#a9d4e7",
    "Lebanon"        = "#ffff33",
    "Iraq"           = "#b15928")

  city_cols <- c(
    "Jeddah"        = "#33a02c",
    "Riyadh"        = "#1f78b4",
    "Dammam"        = "#e31a1c",
    "Qatif"         = "#ff7f00",
    "Muhayil Assir" = "#e7298a",
    "Aljouf"        = "#999999",
    "Sharurah"      = "#5b3a29",
    "Najran"        = "#6a3d9a",
    "Tabouk"        = "#17becf")

  source_shape <- c(Human = "circle", Food = "square", Environment = "triangle")

  mytriangle <- function(coords, v = NULL, params) {
    vcol <- params("vertex", "color")
    if (length(vcol) != 1 && !is.null(v)) vcol <- vcol[v]
    vsz <- params("vertex", "size")
    if (length(vsz) != 1 && !is.null(v)) vsz <- vsz[v]
    r   <- 1.5 * vsz
    ang <- pi/2 + c(0, 2, 4) * pi/3
    for (i in seq_len(nrow(coords))) {
      col_i <- if (length(vcol) == 1) vcol else vcol[i]
      polygon(coords[i, 1] + r * cos(ang), coords[i, 2] + r * sin(ang),
              col = col_i, border = "black", lwd = 0.5)
    }
  }
  add_shape("triangle", clip = shapes("circle")$clip, plot = mytriangle)

  clean      <- function(x) gsub("[^A-Za-z]", "", x)
  loc_lookup <- setNames(location_cols, clean(names(location_cols)))

  snp_width <- function(snp) ifelse(snp <= 1, 6,
                             ifelse(snp <= 5, 4,
                             ifelse(snp <= 10, 2.5, 1.2)))
  sizes     <- list(baps3 = 11)
  node_size <- function(nm) if (!is.null(sizes[[nm]])) sizes[[nm]] else 12
  star_pts <- function(x, y, col, outer = 0.045, inner = 0.018) {
    ang <- pi/2 + seq(0, 2*pi, length.out = 11)[-11]
    rad <- rep(c(outer, inner), 5)
    for (i in seq_along(x)) {
      if (is.na(col[i])) next
      polygon(x[i] + rad * cos(ang), y[i] + rad * sin(ang),
              col = col[i], border = "black", lwd = 0.4)
    }
  }

  make_net <- function(stub, thresh = 20) {
    dna     <- read.dna(sprintf("Figure1/%s.filtered_polymorphic_sites.fasta", stub),
                        format = "fasta")
    distdna <- dist.dna(dna, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
    mf      <- all_meta[match(colnames(distdna), all_meta$label), ]

    ob <- seqTrack(distdna, x.names = colnames(distdna), x.dates = mf$collection_date)
    ob <- ob %>%
      mutate(ances = if_else(weight > thresh & !is.na(weight), NA_real_, ances)) %>%
      drop_na()

    addmat <- matrix(0, nrow(distdna), nrow(distdna),
                     dimnames = list(colnames(distdna), colnames(distdna)))
    snpmat <- addmat
    for (i in seq_len(nrow(ob))) {
      if (!is.na(ob$ances[i])) {
        addmat[ob$ances[i], ob$id[i]] <- 1
        snpmat[ob$ances[i], ob$id[i]] <- ob$weight[i]
      }
    }

    keep <- (colSums(addmat) + rowSums(addmat)) != 0
    net  <- graph_from_adjacency_matrix(addmat[keep, keep, drop = FALSE],
                                        mode = "directed", diag = FALSE)
    el   <- as_edgelist(net)
    E(net)$snp <- snpmat[cbind(el[, 1], el[, 2])]

    mo <- mf[match(V(net)$name, mf$label), ]
    V(net)$location <- mo$location
    V(net)$color    <- loc_lookup[clean(mo$location)]
    V(net)$color[is.na(V(net)$color)] <- "#BDBDBD"
    shp <- source_shape[mo$Source]; shp[is.na(shp)] <- "circle"
    V(net)$shape <- shp
    V(net)$city  <- mo$City
    net
  }

  stubs <- c("baps2", "baps3", "baps4", "baps6", "baps7", "baps8")
  nets  <- lapply(stubs, make_net)
  names(nets) <- stubs

  layout_algo <- list(baps3 = layout_with_kk)
  algo_of     <- function(nm) if (!is.null(layout_algo[[nm]])) layout_algo[[nm]] else layout_with_fr

  nets_prev <- lapply(stubs, make_net, thresh = 10); names(nets_prev) <- stubs
  prev_pos  <- list()
  for (nm in stubs) {
    set.seed(1); pp <- algo_of(nm)(nets_prev[[nm]])
    rownames(pp) <- V(nets_prev[[nm]])$name
    prev_pos[[nm]] <- pp
  }

  gopt <- function(seed) function(g) { set.seed(seed); layout_with_graphopt(g) }
  fresh_layout <- list(baps2 = gopt(1), baps3 = gopt(2), baps7 = gopt(1))

  save_layout <- function(nm, seed = 1) {
    g <- nets[[nm]]
    if (!is.null(fresh_layout[[nm]])) {
      layouts[[nm]] <<- fresh_layout[[nm]](g); return(invisible())
    }
    base <- prev_pos[[nm]]
    m    <- match(V(g)$name, rownames(base)); have <- !is.na(m)
    pos  <- matrix(0, vcount(g), 2); pos[have, ] <- base[m[have], ]
    if (all(have)) { layouts[[nm]] <<- pos; return(invisible()) }
    A <- as_adjacency_matrix(as.undirected(g, mode = "collapse"), sparse = FALSE)
    set.seed(seed)
    for (i in which(!have)) {
      nb <- which(A[i, ] > 0 & have)
      pos[i, ] <- if (length(nb)) colMeans(base[m[nb], , drop = FALSE]) + runif(2, -.1, .1)
                  else runif(2, -.5, .5)
    }
    set.seed(seed)
    layouts[[nm]] <<- layout_with_fr(g, coords = pos, niter = 120, start.temp = 0.10)
  }

  layouts <- list()
  for (nm in stubs) save_layout(nm)

  # Figure 4 - COMBINE: framed panels + Saudi map + text legend
  panel_plot <- function(nm) {
    net <- nets[[nm]]; lay <- layouts[[nm]]; sz <- node_size(nm)
    par(mar = c(0, 0, 2.5, 0))
    pad <- 1.0
    plot(net, layout = lay, rescale = TRUE, xlim = c(-pad, pad), ylim = c(-pad, pad),
         vertex.color = V(net)$color, vertex.shape = V(net)$shape,
         vertex.size  = sz, vertex.label = NA,
         edge.color = "black", edge.arrow.size = 0.5, edge.curved = 0.2,
         edge.width = snp_width(E(net)$snp))
    title(main = toupper(nm), cex.main = 2, font.main = 2)
    rc <- norm_coords(lay, -1, 1, -1, 1)
    star_pts(rc[, 1], rc[, 2], col = city_cols[V(net)$city],
             outer = 0.045 * sz/12, inner = 0.018 * sz/12)

    x <- grconvertX(c(0.006, 0.994), "nfc", "user")
    y <- grconvertY(c(0.006, 0.994), "nfc", "user")
    rect(x[1], y[1], x[2], y[2], border = "black", lwd = 1.5, xpd = NA)
  }

  saudi <- ne_states(country = "Saudi Arabia", returnclass = "sf")
  city_pts <- data.frame(
    city = c("Aljouf", "Qatif", "Dammam", "Riyadh", "Jeddah",
             "Muhayil Assir", "Sharurah", "Najran", "Tabouk"),
    lon  = c(40.2064, 49.9964, 50.0552, 46.6753, 39.1728, 42.0393, 47.1167, 44.1277,
             36.5662),
    lat  = c(29.9697, 26.5652,  26.4257, 24.7136, 21.5433, 18.5444, 17.4667, 17.4924,
             28.3835),
    pos  = c(3, 4, 4, 4, 2, 2, 1, 4, 2),
    stringsAsFactors = FALSE)
  city_pts$col <- city_cols[city_pts$city]

  star_geo <- function(lon, lat, col, R = 0.9) {
    ang <- pi/2 + seq(0, 2*pi, length.out = 11)[-11]
    rad <- rep(c(R, R*0.4), 5); asp <- 1/cos(mean(lat)*pi/180)
    for (i in seq_along(lon))
      polygon(lon[i] + rad*cos(ang)*asp, lat[i] + rad*sin(ang),
              col = col[i], border = "black", lwd = 0.5)
  }

  city_pts$dx <- 0; city_pts$dy <- 0
  city_pts$dx[city_pts$city == "Qatif"]   <- -0.1; city_pts$dy[city_pts$city == "Qatif"]   <-
    0.1
  city_pts$dx[city_pts$city == "Dammam"]  <-  0.1; city_pts$dy[city_pts$city == "Dammam"]  <-
    -0.1
  city_pts$pos[city_pts$city == "Najran"]   <- 1
  city_pts$dy[city_pts$city == "Najran"]    <- -0.8
  city_pts$pos[city_pts$city == "Sharurah"] <- 4
  city_pts$pos[city_pts$city == "Tabouk"] <- 3
  city_pts$pos[city_pts$city == "Riyadh"] <- 2

  city_pts$plon <- city_pts$lon
  city_pts$plat <- city_pts$lat
  city_pts$plon[city_pts$city == "Qatif"]  <- city_pts$plon[city_pts$city == "Qatif"]  - 0.15
  city_pts$plat[city_pts$city == "Qatif"]  <- city_pts$plat[city_pts$city == "Qatif"]  + 1.2
  city_pts$plon[city_pts$city == "Dammam"] <- city_pts$plon[city_pts$city == "Dammam"] + 0.25
  city_pts$plat[city_pts$city == "Dammam"] <- city_pts$plat[city_pts$city == "Dammam"] - 1.3

  map_legend <- function() {
    par(mar = c(10, 0, 0, 9))
    bb   <- st_bbox(saudi)
    padx <- (bb["xmax"] - bb["xmin"]) * 0.10
    pady <- (bb["ymax"] - bb["ymin"]) * 0.10
    plot(st_geometry(saudi), col = "white", border = "grey45", lwd = 0.6,
         xlim = c(bb["xmin"] - padx*2, bb["xmax"] + padx),
         ylim = c(bb["ymin"] - pady, bb["ymax"] + pady))
    star_geo(city_pts$plon, city_pts$plat, city_pts$col, R = 1.75)
    lx <- city_pts$plon + city_pts$dx
    ly <- city_pts$plat + city_pts$dy
    moved <- city_pts$dx != 0 | city_pts$dy != 0
    segments(city_pts$plon[moved], city_pts$plat[moved], lx[moved], ly[moved],
             col = "grey60", lwd = 0.5, xpd = NA)
    text(lx, ly, city_pts$city, pos = city_pts$pos, cex = 1.5, font = 2,
         offset = 0.5, xpd = NA)
  }

  LEFT <- 0.1
  legend_plot <- function() {
    loc <- names(location_cols)[names(location_cols) %in%
                                  unlist(lapply(nets, function(n) V(n)$location))]
    par(mar = c(0, 0, 0, 0)); plot.new(); plot.window(c(0, 1), c(0, 1))
    legend(0, 1.00, loc, pt.bg = location_cols[loc], pch = 21,
           pt.cex = 1.7, cex = 1.5, bty = "n", y.intersp = 0.95,
           title = expression(bold("Country")), title.adj = 0)
    legend(0, 0.52, names(source_shape), pch = c(21, 22, 24), pt.bg = "grey70",
           pt.cex = 1.7, cex = 1.5, bty = "n", y.intersp = 0.95,
           title = expression(bold("Source")), title.adj = 0)
    legend(0, 0.28, c("0-1 SNPs", "2-5 SNPs", "6-10 SNPs", "11-20 SNPs"),
           lwd = c(6, 4, 2.5, 1.2), cex = 1.5, bty = "n", y.intersp = 0.95,
           title = expression(bold("Key")), title.adj = 0)
  }

  assemble_figure4 <- function(file = "Figure1/Figure_4.svg", ncol = 3,
                               width = 18, height = 9, leg_width = 0.8,
                               legend_frac = 0.66, res = 150) {
    use <- names(layouts); n <- length(use); nr <- ceiling(n / ncol)
    K <- 6; R <- nr * K
    m <- matrix(0, R, ncol + 1); idx <- 1
    for (r in seq_len(nr)) for (cc in seq_len(ncol))
      if (idx <= n) { m[((r-1)*K + 1):(r*K), cc] <- idx; idx <- idx + 1 }
    legrows <- max(1, round(R * legend_frac))
    m[1:legrows, ncol + 1]        <- n + 1
    m[(legrows + 1):R, ncol + 1]  <- n + 2

    if (grepl("\\.png$", file, ignore.case = TRUE))
      png(file, width = width * res, height = height * res, res = res)
    else
      svglite::svglite(file, width = width, height = height)

    layout(m, widths = c(rep(1, ncol), leg_width))
    for (nm in use) panel_plot(nm)
    legend_plot()
    map_legend()
    dev.off(); layout(1)
    message("Wrote ", file)
  }

  assemble_figure4("Figure1/Figure_4.svg")

  #Figure5

  library(ape); library(dplyr); library(tidyr); library(readr)
  library(ggtree); library(phytools); library(ggplot2); library(aplot)
  library(RColorBrewer);library(ggtext)

  region_lv  <- c("Middle East", "Europe", "North America", "East Asia",
                  "Oceania", "Africa")
  country_lv <- c("Saudi Arabia", "Kuwait",
                  "United Kingdom",
                  "United States",
                  "China",
                  "Other")
  col_region  <- c("Africa"="#1f78b4","East Asia"="#33a02c","Europe"="#6a3d9a",
                   "Middle East"="#1b9e77","North America"="#ff7f00","Oceania"="#a6761d")
  col_region  <- col_region[region_lv]
  col_source  <- c("Human"="blue","Food"="red","Environment"="#FFD700","NA"="grey70")
  col_country <- c(
    "China"          = "#e78ac3",
    "Kuwait"         = "#a6d854",
    "Saudi Arabia"   = "#8da0cb",
    "United Kingdom" = "#fc8d62",
    "United States"  = "#66c2a5",
    "Other"          = "#b3b3b3")
  col_country <- col_country[country_lv]

  custom_round <- function(x, digits) {
    threshold <- 10^(-digits); rounded <- round(x, digits)
    if (x < threshold & x > 0) rounded <- threshold
    rounded
  }

  top_title_bold <- function(t) list(ggtitle(t),
                                     theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")))
  big_key <- guides(fill = guide_legend(override.aes = list(size = 8)))

  bar_theme <- theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title =
            element_blank(),
          legend.position = "right",
          legend.title = element_text(size = 13, face = "bold"),
          legend.text  = element_text(size = 12),
          legend.key.size = unit(0.35, "cm"),
          legend.spacing.y = unit(0.2, "cm"),
          plot.margin = margin(8, 2, 8, 0))

  dot_theme <- theme_bw() +
    theme(axis.text.x  = element_text(size = 10, angle = 45, vjust = 1, hjust = 1, face =
                                        "bold"),
          axis.text.y  = element_blank(), axis.title.y = element_blank(),
          axis.title.x = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          plot.margin  = margin(10, 25, 10, 30))

  dna     <- read.dna("Figure1/Enteritidis_snp_sites.aln", format = "fasta")
  distdna <- dist.dna(dna, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
  rownames(distdna) <- as.character(rownames(distdna)); colnames(distdna) <-
    rownames(distdna)

  meta <- read_tsv("Figure1/metadata_baps.tsv", show_col_types = FALSE) %>%
    mutate(BAPS = as.character(baps_level2),
           Country = factor(ifelse(location %in% names(col_country), location, "Other"),
                            levels = names(col_country)),
           Source  = ifelse(is.na(Source) | Source == "NA", "NA", Source),
           Source  = factor(Source, levels = c("Human", "Food", "Environment", "NA")),
           region  = factor(region, levels = region_lv))

  clusters <- as.character(sort(unique(as.numeric(meta$BAPS))))
  avg_d <- function(a, b) {
    s1 <- meta$Isolate[meta$BAPS == a]; s2 <- meta$Isolate[meta$BAPS == b]
    mean(distdna[as.character(s1), as.character(s2)], na.rm = TRUE)
  }
  M <- matrix(0, length(clusters), length(clusters), dimnames = list(clusters, clusters))
  for (i in clusters) for (j in clusters) if (i != j) M[i, j] <- avg_d(i, j)
  M <- apply(M, c(1, 2), round)

  njtree_baps <- midpoint.root(bionj(M))
  p1 <- ggtree(njtree_baps) + geom_treescale() +
    geom_tiplab(aes(label = label), align = TRUE, size = 5, offset = 0.2, linesize = .2) +
    xlim(0, max(node.depth.edgelength(njtree_baps)) * 1.1) +
    theme(plot.margin = margin(8, 0, 8, 8))

  lv <- njtree_baps$tip.label
  ny <- length(lv)

  panel_title <- function(lab, xc = 0.5)
    list(annotate("text", x = xc, y = ny + 0.85, label = lab,
                  angle = 0, hjust = 0.5, vjust = 0, size = 4, fontface = "bold"),
         scale_y_discrete(expand = expansion(add = c(0, 0))),
         coord_cartesian(ylim = c(0.5, ny + 0.5), clip = "off"),
         theme(plot.margin = margin(t = 100, r = 2, b = 8, l = 0)))

  pRegion  <- ggplot(meta, aes(1, BAPS, fill = region)) +
    geom_col(position = position_fill(reverse = TRUE), width = 1) +
    scale_fill_manual(name = "Region", values = col_region) + bar_theme + big_key +
    panel_title("Region")

  pCountry <- ggplot(meta, aes(1, BAPS, fill = Country)) +
    geom_col(position = position_fill(reverse = TRUE), width = 1) +
    scale_fill_manual(name = "Country", values = col_country) + bar_theme + big_key +
    panel_title("Country")

  pSource  <- ggplot(meta, aes(1, BAPS, fill = Source)) +
    geom_col(position = position_fill(reverse = TRUE), width = 1) +
    scale_fill_manual(name = "Source", values = col_source) + bar_theme + big_key +
    panel_title("Source")
  all_meta <- read_tsv("Figure1/metadata_baps.tsv", show_col_types = FALSE)
  all_meta$label <- all_meta$Isolate

  amrfinder_result <- read_tsv("Figure1/amrfinder_all.tsv", show_col_types = FALSE) %>%
    filter(`Gene symbol` != "Gene symbol") %>%
    rename(Gene_symbol = `Gene symbol`, Element_type = `Element type`)

  baps_count <- all_meta %>% group_by(BAPS) %>% summarise(sample_amount = n(), .groups = "drop")

  amr_freq <- left_join(all_meta, amrfinder_result, by = c("label" = "Sample_name")) %>%
    filter(Element_type %in% c("AMR")) %>%
    select(label, Gene_symbol, BAPS) %>%
    group_by(BAPS, Gene_symbol) %>% summarise(count = n(), .groups = "drop")

  amr_freq_2 <- amr_freq %>%
    pivot_wider(names_from = Gene_symbol, values_from = count, values_fill = 0) %>%
    pivot_longer(!BAPS, names_to = "gene", values_to = "count") %>%
    left_join(baps_count, by = "BAPS") %>%
    mutate(frequency = count / sample_amount) %>%
    filter(frequency != 0, gene != "blaTEM")

  gene_levels <- c(
    "aac(3)-IId","aadA22","aph(3'')-Ib","aph(3')-IIa","aph(6)-Id",
    "blaCTX-M-8","blaSHV-12","blaTEM-1",
    "acrB_R717L","mdsA","mdsB","ramR_T18P",
    "gyrA_D87Y","gyrA_S83Y","qnrB19","qnrS13",
    "dfrA14","sul2","mph(A)","floR",
    "tet(A)","lnu(F)","ble")
  italicize_gene <- function(x) {
    ifelse(grepl("_", x),
           paste0("*", sub("_.*", "", x), "* ", sub(".*_", "", x)),
           paste0("*", x, "*"))
  }
  gene_labs <- setNames(italicize_gene(gene_levels), gene_levels)

  drug_class_map <- c(
    "aac(3)-IId"="Aminoglycosides","aadA22"="Aminoglycosides","aph(3'')-Ib"="Aminoglycosides",
    "aph(3')-IIa"="Aminoglycosides","aph(6)-Id"="Aminoglycosides",
    "blaCTX-M-8"="Beta-lactams","blaSHV-12"="Beta-lactams","blaTEM-1"="Beta-lactams",
    "acrB_R717L"="Multidrug efflux","mdsA"="Multidrug efflux","mdsB"="Multidrug efflux",
    "ramR_T18P"="Multidrug efflux",
    "gyrA_D87Y"="Fluoroquinolones","gyrA_S83Y"="Fluoroquinolones",
    "qnrB19"="Fluoroquinolones","qnrS13"="Fluoroquinolones",
    "dfrA14"="Trimethoprim","sul2"="Sulfonamides","mph(A)"="Macrolides",
    "floR"="Phenicols","tet(A)"="Tetracyclines","lnu(F)"="Lincosamides",
    "ble"="Bleomycin")

  class_levels <- unique(drug_class_map[gene_levels])

  amr_freq_2$BAPS <- factor(as.character(amr_freq_2$BAPS), levels = lv)
  amr_freq_2$gene <- factor(amr_freq_2$gene, levels = gene_levels)

  amr_freq_2$drug_class <- factor(drug_class_map[as.character(amr_freq_2$gene)], levels = class_levels)

  tot_n    <- sum(baps_count$sample_amount)
  gene_tot <- amr_freq %>% group_by(Gene_symbol) %>%
    summarise(gene_count = sum(count), .groups = "drop")

  amr_sig <- amr_freq_2 %>%
    left_join(gene_tot, by = c("gene" = "Gene_symbol")) %>%
    rowwise() %>%
    mutate(
      p_value = fisher.test(matrix(c(count, sample_amount - count,
                                     gene_count - count,
                                     (tot_n - sample_amount) - (gene_count - count)),
                                   nrow = 2, byrow = TRUE),
                            alternative = "greater")$p.value) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           sig   = !is.na(p_adj) & p_adj < 0.05)

  amr_sig$drug_class <- factor(drug_class_map[as.character(amr_sig$gene)], levels = class_levels)

  class_annot <- tibble::tibble(gene = levels(droplevels(amr_freq_2$gene))) %>%
    mutate(x = dplyr::row_number(),
           drug_class = drug_class_map[gene]) %>%
    group_by(drug_class) %>%
    summarise(xmin = min(x), xmax = max(x), .groups = "drop")

  ny <- length(unique(amr_freq_2$BAPS))

  gene_location <- c(
    "acrB_R717L"="Chromosome","mdsA"="Chromosome","mdsB"="Chromosome","ramR_T18P"="Chromosome",
    "gyrA_D87Y"="Chromosome","gyrA_S83Y"="Chromosome",
    "blaTEM-1"="Plasmid","blaCTX-M-8"="Plasmid","blaSHV-12"="Plasmid",
    "qnrB19"="Plasmid","qnrS13"="Plasmid","aadA22"="Plasmid","lnu(F)"="Plasmid",
    "tet(A)"="Plasmid","floR"="Plasmid","aph(6)-Id"="Plasmid","aph(3'')-Ib"="Plasmid",
    "aac(3)-IId"="Unknown","aph(3')-IIa"="Unknown","dfrA14"="Unknown",
    "sul2"="Unknown","mph(A)"="Unknown","ble"="Unknown")
  gene_loc_df <- tibble::tibble(gene = levels(droplevels(amr_freq_2$gene))) %>%
    dplyr::mutate(location = factor(gene_location[gene],
                                    levels = c("Chromosome","Plasmid","Unknown")))

  y_loc  <- ny + 0.75
  yline  <- ny + 1.02
  ylabel <- ny + 1.12

  p2 <- ggplot(amr_freq_2, aes(gene, BAPS)) +
    geom_tile(data = dplyr::filter(amr_sig, sig), aes(gene, BAPS),
              inherit.aes = FALSE, fill = "red", alpha = 0.22,
              width = 0.9, height = 0.9) +

    geom_point(aes(size = frequency), alpha = 0.7,
               color = "#b4a6d4", fill = "#b4a6d4", shape = 21) +

    geom_text(aes(label = sapply(frequency, custom_round, digits = 2)),
              size = 3.5) +

    geom_segment(data = class_annot, inherit.aes = FALSE,
                 aes(x = xmin - 0.4, xend = xmax + 0.4,
                     y = yline, yend = yline),
                 linewidth = 0.5, colour = "grey30") +

    geom_text(data = class_annot, inherit.aes = FALSE,
              aes(x = xmin - 0.4, y = ylabel, label = drug_class),
              angle = 45, hjust = 0, vjust = 0, size = 9 / .pt) +

    geom_tile(data = gene_loc_df, inherit.aes = FALSE,
              aes(x = gene, y = y_loc, fill = location),
              width = 0.9, height = 0.4, colour = "white",, linewidth = 0.25) +
    scale_fill_manual(values = c(Chromosome = "#8FA6BF", Plasmid = "#D5A07A", Unknown ="#D2D2D2"),
                      name = "Gene/Mutation Location", drop = FALSE,
                      guide = guide_legend(title.theme = element_text(face = "bold", size = 13))) +

    scale_size(range = c(0.1, 8), limits = c(0, 1),
               breaks = c(0.25, 0.5, 0.75, 1),
               name = "Frequency",
               guide = guide_legend(
                 title.theme = element_text(face = "bold", size = 13))) +
    scale_x_discrete(labels = gene_labs) +

    scale_y_discrete(expand = expansion(add = c(0, 0))) +

    coord_cartesian(ylim = c(0.5, ny + 0.5), clip = "off") +

    labs(x = "Antibiotic Resistance Determinants") +

    dot_theme +
    theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, size = 9),
          axis.title.x = element_text(size = 13, face = "bold"),
          plot.margin = margin(t = 100, r = 40, b = 10, l = 10),
          legend.title    = element_text(size = 13, face = "bold"),
          legend.text     = element_text(size = 12),
          legend.justification = "top")

  baps_plasmid <- read.csv("Figure1/plasmid_freq.csv", header = TRUE) %>%
    mutate(across(where(is.numeric), ~ if_else(. < 90, 0L, 1L)))
  baps_plasmid <- left_join(all_meta, baps_plasmid, by = c("label" = "Sample"))
  for (cl in c("SE016_3","SE055_2","SE055_3","SE061_2","SE119_2","SE242_2"))
    baps_plasmid[[cl]] <- as.numeric(baps_plasmid[[cl]])
  baps_plasmid$BAPS <- as.numeric(baps_plasmid$BAPS)

  plasmid_freq <- baps_plasmid %>% group_by(BAPS) %>%
    summarise(SE016_3 = sum(SE016_3, na.rm = TRUE),
              SE061_2 = sum(SE061_2, na.rm = TRUE),
              SE055_3 = sum(SE055_3, na.rm = TRUE),
              SE055_2 = sum(SE055_2, na.rm = TRUE),
              SE119_2 = sum(SE119_2, na.rm = TRUE),
              SE242_2 = sum(SE242_2, na.rm = TRUE),
              count = n(), .groups = "drop") %>%
    mutate(across(-c(BAPS, count), ~ . / count)) %>%
    select(-count) %>%
    pivot_longer(!BAPS, names_to = "plasmid", values_to = "frequency")

  plasmid_freq$BAPS <- factor(as.character(plasmid_freq$BAPS), levels = lv)
  plasmid_freq$plasmid <- factor(plasmid_freq$plasmid)

  .rep <- read.csv("Figure1/plasmid_replicon_by_contig.csv")
  replicon_lu <- setNames(.rep$replicon_type, .rep$contig)
  replicon_annot <- data.frame(
    x        = seq_along(levels(plasmid_freq$plasmid)),
    replicon = replicon_lu[levels(plasmid_freq$plasmid)])

  ny  <- length(lv)
  y_plasmid <- ny + 2.0
  p4 <- ggplot(plasmid_freq, aes(plasmid, BAPS, size = frequency)) +
    geom_point(data = subset(plasmid_freq, frequency > 0),
               alpha = 0.7,
               color = "#b4a6d4",
               fill = "#b4a6d4",
               shape = 21) +
    geom_text(data = subset(plasmid_freq, frequency > 0),
              aes(label = sapply(frequency, custom_round, digits = 2)),
              size = 3.5) +
    geom_text(data = replicon_annot, inherit.aes = FALSE,
              aes(x = x, y = ny + 0.6, label = replicon),
              angle = 45, hjust = 0, vjust = 0, size = 9 / .pt) +
    scale_size(range = c(0.1, 8),
               limits = c(0, 1),
               guide = "none") +
    scale_y_discrete(expand = expansion(add = c(0, 0))) +
    coord_cartesian(ylim = c(0.5, ny + 0.5), clip = "off") +
    labs(x = "Plasmids") +
    dot_theme +
    theme(
      axis.title.x = element_text(size = 13, face = "bold"),
      plot.margin = margin(t = 150, r = 40, b = 10, l = 10)
    )
  prophage_cov <- read_tsv("Figure1/prophage_freq.tsv", show_col_types = FALSE)
  prophage_cov[is.na(prophage_cov)] <- 0
  baps_prophage <- left_join(all_meta, prophage_cov, by = c("label" = "Sample"))
  baps_prophage$`ST64B (AY055382.1)` <- as.numeric(baps_prophage$`ST64B (AY055382.1)`)
  baps_prophage$BAPS <- as.numeric(baps_prophage$BAPS)

  prophage_freq <- baps_prophage %>% group_by(BAPS) %>%
    summarise(ST64B = sum(`ST64B (AY055382.1)`, na.rm = TRUE), count = n(), .groups = "drop") %>%
    mutate(ST64B = ST64B / count) %>% select(-count) %>%
    pivot_longer(!BAPS, names_to = "prophage", values_to = "frequency")
  prophage_freq$BAPS <- factor(as.character(prophage_freq$BAPS), levels = lv)
  prophage_freq$prophage <- factor(prophage_freq$prophage)

  ny  <- length(lv)
  y_prophage <- ny + 0.85
  p6 <- ggplot(prophage_freq, aes(prophage, BAPS, size = frequency)) +
    geom_point(alpha = 0.7,
               color = "#b4a6d4",
               fill = "#b4a6d4",
               shape = 21) +
    geom_text(aes(label = sapply(frequency, custom_round, digits = 2)),
              size = 3.5) +
    scale_size(range = c(0.1, 8),
               limits = c(0, 1),
               guide = "none") +
    scale_y_discrete(
      limits = lv,
      drop = FALSE,
      expand = expansion(add = c(0, 0))
    ) +
    coord_cartesian(ylim = c(0.5, ny + 0.5), clip = "off") +
    labs(x = "Prophage") +
    dot_theme +
    theme(
      axis.title.x = element_text(size = 13, face = "bold"),
      plot.margin = margin(t = 100, r = 10, b = 10, l = 10)
    )

  library(cowplot)

  top_space <- ggplot() + theme_void()

  extract_legend <- function(p) {
    g   <- ggplot2::ggplotGrob(p)
    nms <- vapply(g$grobs, function(x) if (is.null(x$name)) "" else x$name, character(1))
    cand <- g$grobs[grepl("guide-box", nms)]
    for (b in cand) {
      if (!inherits(b, "zeroGrob") && !is.null(b$grobs) && length(b$grobs) > 0) return(b)
    }
    if (length(cand)) cand[[1]] else ggplot2::zeroGrob()
  }

  legs <- list(
    extract_legend(pRegion),
    extract_legend(pCountry),
    extract_legend(pSource),
    extract_legend(p2)
  )
  wds  <- vapply(legs, function(g) sum(grid::convertWidth(grid::grobWidth(g), "cm", TRUE)), numeric(1))
  legs <- Map(function(g, w)
    if (max(wds) - w > 0.01) gtable::gtable_add_cols(g, grid::unit(max(wds) - w, "cm")) else g,
    legs, wds)
  leg_h <- vapply(legs, function(g)
    sum(grid::convertHeight(grid::grobHeight(g), "cm", valueOnly = TRUE)), numeric(1))

  noleg <- theme(legend.position = "none")
  combined <- (pRegion + noleg) %>%
    insert_left(p1,                width = 2.5) %>%
    insert_right(pCountry + noleg, width = 1) %>%
    insert_right(pSource  + noleg, width = 1) %>%
    insert_right(p2       + noleg, width = 8) %>%
    insert_right(p4       + noleg, width = 2.3) %>%
    insert_right(p6       + noleg, width = 0.5) %>%
    insert_top(top_space,          height = 0.25)

  combined_grob <- grid::grid.grabExpr(print(combined))

  FIG_H   <- 9
  TOP_GAP <- 4
  lh      <- leg_h * 1.1
  BOTTOM  <- max(0.5, FIG_H * 2.54 - TOP_GAP - sum(lh))
  legend_block <- cowplot::plot_grid(
    plotlist    = c(list(NULL), legs, list(NULL)),
    ncol        = 1,
    rel_heights = c(TOP_GAP, lh, BOTTOM))

  final <- cowplot::plot_grid(combined_grob, NULL, legend_block,
                              ncol = 3, rel_widths = c(1, 0.015, 0.16))
  final

  ggsave("Figure1/Figure_5.svg", final, width = 15, height = 9,
         device = svglite::svglite)
