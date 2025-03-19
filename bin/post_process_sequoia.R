suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(data.table))
suppressMessages(library(ggtree))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(ape))
suppressMessages(library(paletteer))

###### Input files #######
###### Here substitute for the files located in your Sequoia output folder #####
args <- list(
  sequoia_dir = commandArgs(trailingOnly = TRUE)[1],
  out_dir = commandArgs(trailingOnly = TRUE)[2])

file_placed_variants <- list.files(args$sequoia_dir, pattern = "assigned_to_branches.txt", full.names = TRUE)
file_nv <- list.files(args$sequoia_dir, pattern = "NV_filtered_all.txt", full.names = TRUE)
file_nr <- list.files(args$sequoia_dir, pattern = "NR_filtered_all.txt", full.names = TRUE)
file_tree <- list.files(args$sequoia_dir, pattern = "tree_with_branch_length.tree", full.names = TRUE)

theme_ohchibi <- function(size_axis_text.x = 12, size_axis_text.y = 12, size_axis_title.x = 13,
                          size_axis_title.y = 13, angle_text.x = 0, angle_text.y = 0,
                          legend_proportion_size = 1, size_title_text = 13, size_legend_text = 12,
                          size_lines_panel = 0.3, size_panel_border = 0.75, x_hjust = 0.5,
                          y_vjust = 0.5, font_family = "Helvetica", font_face = "plain",
                          size_ticks = 0.5) {
  theme(
    strip.background.x = element_blank(), strip.background.y = element_blank(),
    strip.text.x = element_text(size = size_axis_title.x),
    strip.text.y = element_text(size = size_axis_title.y),
    panel.background = element_rect(fill = "white"), panel.grid.major.x = element_line(
      color = "grey89",
      size = 0.25
    ), panel.grid.major.y = element_line(
      color = "grey89",
      size = 0.25
    ), panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(), panel.border = element_rect(
      fill = NA,
      color = "black", size = size_panel_border
    ), axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", size = size_ticks),
    axis.text.x = element_text(
      family = font_family, face = font_face,
      size = size_axis_text.x, colour = "black", hjust = x_hjust,
      angle = angle_text.x
    ), axis.text.y = element_text(
      family = font_family,
      face = font_face, size = size_axis_text.y, colour = "black",
      vjust = y_vjust, angle = angle_text.y
    ), axis.title.x = element_text(
      family = font_family,
      face = font_face, size = size_axis_title.x, colour = "black"
    ),
    axis.title.y = element_text(
      family = font_family, face = font_face,
      size = size_axis_title.y, colour = "black"
    ), legend.background = element_blank(),
    legend.key.size = unit(legend_proportion_size, "line"),
    legend.title = element_text(
      size = size_title_text, family = font_family,
      face = font_face, colour = "black"
    ), legend.key = element_blank(),
    legend.text = element_text(
      size = size_legend_text, family = font_family,
      face = font_face, colour = "black"
    ), legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = size_axis_title.x)
  )
}

oh.save.pdf <- function(p = NULL, outname = NULL, outdir = "figures/", width = 20,
                        height = 15, family = "Arial", fallback_resolution = 1200,
                        antialias = "default", pointsize = 12) {
  dir.create(outdir, showWarnings = FALSE)
  myfilename <- paste(outdir, outname, sep = "/")
  cairo_pdf(
    filename = myfilename, onefile = FALSE, fallback_resolution = fallback_resolution,
    width = width, height = height, family = family, antialias = antialias,
    pointsize = pointsize
  )
  print(p)
  dev.off()
}


### Read information about mutations matching
df_placed <- read.table(file = file_placed_variants, header = TRUE) %>%
  dplyr::mutate(.data = ., VariantId = paste(Chr, Pos, Ref, Alt, sep = "_"))


# Determine the number of variants per branch
df_freq <- df_placed$Branch %>%
  table() %>%
  sort(decreasing = TRUE) %>%
  data.frame() %>%
  dplyr::rename(.data = ., node = ".")


df_freq$Log2Freq <- log2(df_freq$Freq)

df_freq$node <- as.character(df_freq$node) %>% as.numeric()


# Read nv matrix and melt so we can merge and estimate VAF
df_nv <- read.table(file_nv, header = TRUE, check.names = FALSE) %>%
  as.matrix() %>%
  reshape2::melt()

df_nr <- read.table(file_nr, header = TRUE, check.names = FALSE) %>%
  as.matrix() %>%
  reshape2::melt()

df_vaf <- merge(df_nv, df_nr, by = c("Var1", "Var2")) %>%
  dplyr::mutate(.data = ., VAF = value.x / value.y) %>%
  dplyr::rename(.data = ., VariantId = Var1, biosampleName = Var2, NV = value.x, NR = value.y) %>%
  dplyr::filter(.data = ., VariantId %in% df_placed$VariantId)

df_vaf$PlacedNodeId <- match(df_vaf$VariantId, df_placed$VariantId) %>% df_placed$Branch[.]

rm(df_nv)
rm(df_nr)
gc()


## Script identifying the nodes ##
tree <- ape::read.tree(file = file_tree)

tree2 <- ape::ladderize(tree, right = TRUE)
is_internal <- tree2$edge[, 2] > length(tree2$tip.label)
order_nodes <- tree2$edge[is_internal, 2]


df_nodes_ids <- tidytree::as_tibble(ape::makeNodeLabel(tree, method = "md5sum")) %>% as.data.frame()

###### Here we output the table with the translation

write.table(df_nodes_ids, file = file.path(args$out_dir, "df_translation_nodes_to_labels.tsv"), quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

df_vaf$PlacedNodeLabel <- match(df_vaf$PlacedNodeId, df_nodes_ids$node) %>% df_nodes_ids$label[.]


tree$edge.length <- rep(1, nrow(tree$edge))
tree <- full_join(tree, df_freq, by = "node")


p_tree <- ggtree(tr = tree, aes(color = Log2Freq), ladderize = FALSE, size = 2) +
  scale_color_paletteer_c("viridis::plasma",
    na.value = "#D9D9D9",
    oob = squish,
    breaks = c(0, 5, 10, 100),
    labels = c(0, 32, 1000, 20000),
    name = "Number of\nplaced mutations"
  ) +

  geom_tiplab(size = 0, align = TRUE) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 0.5),
    legend.title = element_text(size = 9),
  )

p_tree$data$InternalNode <- p_tree$data$node

p_tree <- p_tree +
  geom_label(aes(label = InternalNode), color = "black", size = 2.5) +
  geom_tiplab(size = 0, align = TRUE, color = "black") +
  coord_cartesian(clip = "off")

# Get the order of the tips
extract_tree_data <- function(tree_disp, displayorder = TRUE) {
  td_out <- tree_disp$data
  if (displayorder) {
    td_out <- dplyr::arrange(td_out, y)
  }
  return(td_out)
}

order_tips_displayed <- extract_tree_data(p_tree) %>%
  subset(isTip == TRUE) %>%
  as.data.frame() %$% InternalNode


order_labels_displayed <- extract_tree_data(p_tree) %>%
  subset(isTip == TRUE) %>%
  as.data.frame() %$% label


### Now create the VAF graph for all of them
mnums <- df_placed[, c("VariantId", "Branch")] %>%
  unique() %$% Branch %>%
  table()

df_numplaced <- data.frame(Branch = names(mnums), NumPlaced = mnums %>% as.numeric())

df_numplaced$Label <- paste0(df_numplaced$Branch, " (n=", df_numplaced$NumPlaced, ")")


df_vaf$NumPlacedVariantsInNodeId <- match(df_vaf$PlacedNodeId, df_numplaced$Branch) %>% df_numplaced$NumPlaced[.]
df_vaf$LabelFacet <- match(df_vaf$PlacedNodeId, df_numplaced$Branch) %>% df_numplaced$Label[.]


df_vaf_int <- df_vaf %>%
  dplyr::filter(.data = ., PlacedNodeId %in% order_nodes) %>%
  droplevels()

df_vaf_int$PlacedNodeId <- df_vaf_int$PlacedNodeId %>% factor(levels = order_nodes)

df_vaf_int$biosampleName <- df_vaf_int$biosampleName %>% factor(levels = order_labels_displayed)

# Add the actual number of variants to the graph
order_nl <- with(df_vaf_int, order(PlacedNodeId)) %>%
  df_vaf_int$LabelFacet[.] %>%
  unique()

df_vaf_int$LabelFacet <- df_vaf_int$LabelFacet %>% factor(levels = order_nl)

# ORder by chromosome
df_vaf_int$VAF[which(df_vaf_int$VAF == 0)] <- NA

p_vaf_branches <- ggplot(data = df_vaf_int, aes(VariantId, biosampleName)) +
  geom_raster(aes(fill = VAF)) +
  facet_grid(. ~ LabelFacet, space = "fixed", scales = "free") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_line(size = 0.2),
    strip.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 0),
    strip.text.y = element_text(size = 8),
    legend.text = element_text(size = 9, angle = 90, vjust = 0, hjust = 0),
    legend.title = element_text(size = 9),
    plot.title = element_text(size = 9),
  ) +
  scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71", na.value = "#F5F5F5") +
  xlab(label = "Somatic placed variants") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

composition <- egg::ggarrange(p_tree, p_vaf_branches, nrow = 1, widths = c(0.4, 1), draw = FALSE)
dev.off()

# Each sample is approximately 0.4 in width
tot_width <- (df_vaf$biosampleName %>% unique() %>% as.character() %>% length()) * 0.4

# Each sample is approximately 0.2111111 in heigth
tot_height <- (df_vaf$biosampleName %>% unique() %>% as.character() %>% length()) * 0.2111111


oh.save.pdf(p = composition, outname = "res_composition.pdf", outdir = args$out_dir, width = tot_width, height = tot_height)
