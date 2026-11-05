# Publication-ready HLA–KIR allele Manhattan plot
# Notebook:
#   %%R -i data -i bonferonni_threshold -i label_phenotype -i ancestry -i h -i k -i t
# Optional extra input: -i phenotype_list (phenotype id -> description).
#
# data columns: MarkerID, phenotype, neglog10p. p_fdr orders labels.
#   p_bonferroni, when present, selects which points are annotated.
# ancestry: cohort label for the title and file name. A length-1 `phenotype`
#   is used only when ancestry is absent (older cells).
# bonferonni_threshold: y position of the dashed line, on the -log10(P) scale.
#   A p-value in (0, 1) is converted with -log10. `bonferroni_threshold` is
#   accepted as an alias. Either may be a list with neglog10_threshold or
#   p_threshold.
# label_phenotype: phenotype ids or descriptions to annotate (at most 10,
#   smallest p_fdr). Omit it to annotate any phenotype.
# t: trait group in the title and file name. h and k are file-name tags only.
#   `t` is captured before use so it cannot fall through to base::t.

library(ggplot2)
library(dplyr)
library(ggrepel)

scalar_chr <- function(name) {
  if (!exists(name, inherits = TRUE)) return("")
  v <- get(name, inherits = TRUE)
  if (is.function(v) || is.null(v) || !is.atomic(v) || length(v) != 1L) return("")
  v <- trimws(as.character(v))
  if (length(v) != 1L || is.na(v)) "" else v
}

as_neglog10_threshold <- function(raw) {
  if (is.null(raw)) return(5)
  if (is.list(raw) || is.data.frame(raw)) {
    nm <- names(raw)
    if ("neglog10_threshold" %in% nm) raw <- raw[["neglog10_threshold"]]
    else if ("p_threshold" %in% nm) raw <- raw[["p_threshold"]]
    else stop("bonferonni_threshold needs neglog10_threshold or p_threshold")
  }
  raw <- suppressWarnings(as.numeric(raw)[1])
  if (length(raw) != 1L || is.na(raw) || raw <= 0) {
    stop("bonferonni_threshold must be a positive number")
  }
  # Values below 1 are p-values. A -log10 line is at least ~1.3 for alpha 0.05.
  if (raw < 1) raw <- -log10(raw)
  raw
}

if (!exists("data") || !is.data.frame(data)) {
  stop("`data` must be a data frame with MarkerID, phenotype, and neglog10p")
}
if (!"neglog10p" %in% names(data)) {
  pcol <- intersect(c("p", "p.value", "p_value", "pval", "pvalue", "P"), names(data))
  if (!length(pcol)) stop("data needs neglog10p or a p-value column")
  data$neglog10p <- -log10(as.numeric(data[[pcol[1]]]))
}
missing_cols <- setdiff(c("MarkerID", "phenotype", "neglog10p"), names(data))
if (length(missing_cols)) {
  stop("data is missing columns: ", paste(missing_cols, collapse = ", "))
}

ANC <- scalar_chr("ancestry")
if (!nzchar(ANC)) ANC <- scalar_chr("phenotype")
trait <- scalar_chr("t")
h_tag <- scalar_chr("h")
k_tag <- scalar_chr("k")

threshold <- if (exists("bonferonni_threshold", inherits = TRUE)) {
  as_neglog10_threshold(get("bonferonni_threshold", inherits = TRUE))
} else if (exists("bonferroni_threshold", inherits = TRUE)) {
  as_neglog10_threshold(get("bonferroni_threshold", inherits = TRUE))
} else {
  5
}

okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

# Classical HLA order, then KIR centromeric -> telomeric (same order as
# hla_linkage_plot.KIR_ORDER). Genes not listed are appended alphabetically.
canonical_gene_order <- c(
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
  "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
  "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
  "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
  "KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DL5B", "KIR2DS3",
  "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR3DS1",
  "KIR2DL5A", "KIR2DS5", "KIR2DS1", "KIR2DS4", "KIR3DL2"
)

data$MarkerID <- trimws(as.character(data$MarkerID))
data$phenotype <- as.character(data$phenotype)
data$gene <- sub("\\*.*$", "", data$MarkerID)
data$gene <- trimws(sub("\\.+$", "", data$gene))

present_genes <- unique(data$gene[nzchar(data$gene)])
ordered_known <- canonical_gene_order[canonical_gene_order %in% present_genes]
ordered_rest <- sort(setdiff(present_genes, canonical_gene_order))
gene_levels <- c(ordered_known, ordered_rest)

# One x position per allele. Rows are allele x phenotype, so position is the
# allele index within the gene, not the row index.
data_sorted <- data %>%
  filter(nzchar(gene), is.finite(neglog10p)) %>%
  mutate(gene = factor(gene, levels = gene_levels)) %>%
  arrange(gene, MarkerID, phenotype) %>%
  group_by(gene) %>%
  mutate(position = match(MarkerID, unique(MarkerID))) %>%
  ungroup()

if (!nrow(data_sorted)) stop("no finite allele associations to plot")

gene_order <- levels(droplevels(data_sorted$gene))
n_genes <- length(gene_order)

# ---- Gene-block layout -------------------------------------------------
# Each gene gets an x-axis slot whose width tracks its allele count, so genes
# with many alleles are not squashed into the same width as single-allele
# genes. Slots also carry a minimum width, derived from the figure width and
# the tick font size, so neighbouring gene tick labels can never collide.

PLOT_WIDTH_IN <- 8
PLOT_HEIGHT_IN <- 4
AXIS_LABEL_PT <- 10

usable_width_pt <- (PLOT_WIDTH_IN - 1.1) * 72

label_footprint_pt <- function(pt) pt * 1.35 + 3

gene_sizes <- data_sorted %>%
  group_by(gene) %>%
  summarise(n_alleles = n_distinct(MarkerID), .groups = "drop") %>%
  mutate(gene = as.character(gene))
sizes <- gene_sizes$n_alleles[match(gene_order, gene_sizes$gene)]
if (anyNA(sizes) || any(sizes < 1L)) stop("could not count alleles for every gene")

axis_label_pt <- AXIS_LABEL_PT
while (axis_label_pt > 5 &&
       label_footprint_pt(axis_label_pt) * n_genes > 0.9 * usable_width_pt) {
  axis_label_pt <- axis_label_pt - 0.5
}

gap <- max(2, 0.08 * max(sizes))

need_frac <- label_footprint_pt(axis_label_pt) / usable_width_pt
min_slot <- 1
for (i in seq_len(200)) {
  total_width <- sum(pmax(sizes, min_slot) + gap)
  target <- need_frac * total_width - gap
  if (target <= min_slot + 1e-6) break
  min_slot <- min_slot + 0.6 * (target - min_slot)
}
min_slot <- ceiling(min_slot)

slot_width <- pmax(sizes, min_slot)
gene_offsets <- data.frame(
  gene = gene_order,
  slot_width = slot_width,
  offset = cumsum(c(0, head(slot_width + gap, -1L))),
  pad = (slot_width - sizes) / 2,
  stringsAsFactors = FALSE
)

# Strongest associations last, so they are drawn on top of the shared x.
plot_data <- data_sorted %>%
  mutate(gene = as.character(gene)) %>%
  left_join(gene_offsets, by = "gene") %>%
  mutate(
    bp_cum = offset + pad + position - 0.5,
    gene = factor(gene, levels = gene_order),
    gene_short = sub("^(HLA|KIR)-", "", as.character(gene))
  ) %>%
  arrange(neglog10p)

axis_set <- gene_offsets %>%
  mutate(
    gene_short = sub("^(HLA|KIR)-", "", gene),
    center = offset + slot_width / 2
  ) %>%
  arrange(offset)

gene_colors <- setNames(rep(okabe_ito, length.out = n_genes), gene_order)

# phenotype id -> display name. PheWAS pheinfo uses `description`.
pheno_key <- NULL
if (exists("phenotype_list", inherits = TRUE)) {
  phenotype_list <- get("phenotype_list", inherits = TRUE)
  if (!is.null(phenotype_list) && !is.atomic(phenotype_list)) {
    pl <- as.data.frame(phenotype_list, stringsAsFactors = FALSE)
    desc_col <- intersect(
      c("description", "phecode_string", "label", "phenotype_description"),
      names(pl)
    )
    if ("phenotype" %in% names(pl) && length(desc_col)) {
      pheno_key <- pl %>%
        transmute(
          phenotype = as.character(phenotype),
          phenotype_label = as.character(.data[[desc_col[1]]])
        ) %>%
        distinct(phenotype, .keep_all = TRUE)
    }
  }
}

wanted <- NULL
if (exists("label_phenotype", inherits = TRUE)) {
  label_phenotype <- get("label_phenotype", inherits = TRUE)
  if (!is.null(label_phenotype)) {
    if (is.data.frame(label_phenotype)) {
      col <- if ("phenotype" %in% names(label_phenotype)) "phenotype" else names(label_phenotype)[1]
      wanted <- as.character(label_phenotype[[col]])
    } else {
      wanted <- as.character(label_phenotype)
    }
    wanted <- wanted[!is.na(wanted) & nzchar(wanted)]
    if (!is.null(pheno_key) && length(wanted)) {
      wanted <- unique(c(
        wanted,
        pheno_key$phenotype[pheno_key$phenotype_label %in% wanted],
        pheno_key$phenotype_label[pheno_key$phenotype %in% wanted]
      ))
    }
    if (!length(wanted)) wanted <- character()
  }
}

label_data <- plot_data
if ("p_bonferroni" %in% names(label_data)) {
  label_data <- label_data %>% filter(is.finite(p_bonferroni), p_bonferroni < 0.05)
} else {
  label_data <- label_data %>% filter(neglog10p >= threshold)
}

if (!is.null(pheno_key)) {
  label_data <- label_data %>%
    left_join(pheno_key, by = "phenotype") %>%
    mutate(
      phenotype_label = if_else(
        is.na(phenotype_label) | !nzchar(phenotype_label),
        phenotype,
        phenotype_label
      )
    )
} else {
  label_data <- label_data %>% mutate(phenotype_label = phenotype)
}

if (!is.null(wanted)) {
  label_data <- label_data %>%
    filter(phenotype %in% wanted | phenotype_label %in% wanted)
}

if ("p_fdr" %in% names(label_data)) {
  label_data <- label_data %>% arrange(p_fdr)
} else {
  label_data <- label_data %>% arrange(desc(neglog10p))
}

label_data <- label_data %>%
  distinct(MarkerID, phenotype, .keep_all = TRUE) %>%
  slice_head(n = 10L) %>%
  mutate(
    # End-anchored novel tag only. A bare gsub("ne") deletes those letters
    # inside the allele name.
    allele_clean = sub("(:?new)$", "", MarkerID),
    allele_clean = trimws(sub("\\.+$", "", allele_clean)),
    label = paste0(allele_clean, "\n", phenotype_label)
  )

title_text <- paste(c(
  if (nzchar(ANC)) ANC,
  "HLA-KIR",
  if (nzchar(trait)) trait,
  "trait association"
), collapse = " ")

manhattan_plot_obj <- ggplot(plot_data, aes(x = bp_cum, y = neglog10p, color = gene)) +
  geom_point(alpha = 0.85, size = 1.8, stroke = 0) +
  geom_hline(yintercept = threshold, color = "#D55E00", linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(
    breaks = axis_set$center,
    labels = axis_set$gene_short,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  scale_color_manual(values = gene_colors, drop = FALSE) +
  labs(
    title = title_text,
    x = NULL,
    y = expression(-log[10](italic(P)))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = axis_label_pt, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    plot.margin = margin(8, 12, 8, 8)
  )

if (nrow(label_data)) {
  manhattan_plot_obj <- manhattan_plot_obj +
    geom_text_repel(
      data = label_data,
      aes(label = label),
      color = "black",
      size = 3.2,
      lineheight = 0.9,
      box.padding = 0.45,
      point.padding = 0.3,
      min.segment.length = 0,
      segment.color = "grey40",
      segment.size = 0.3,
      max.overlaps = Inf,
      seed = 42,
      show.legend = FALSE
    )
}

file_token <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9._+-]+", "_", x))
name_parts <- c(ANC, "HLA", trait, h_tag, k_tag, "manhattan")
name_parts <- file_token(name_parts)
out_base <- paste(name_parts[nzchar(name_parts)], collapse = "_")

ggsave(paste0(out_base, ".pdf"), plot = manhattan_plot_obj,
       width = PLOT_WIDTH_IN, height = PLOT_HEIGHT_IN, units = "in")
ggsave(paste0(out_base, ".png"), plot = manhattan_plot_obj,
       width = PLOT_WIDTH_IN, height = PLOT_HEIGHT_IN, units = "in", dpi = 300)

manhattan_plot_obj
