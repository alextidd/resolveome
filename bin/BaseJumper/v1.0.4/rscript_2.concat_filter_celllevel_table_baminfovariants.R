library(dplyr)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)

cutoff_as <- args[1] %>% as.numeric()
cutoff_clip <- args[2] %>% as.numeric()
cutoff_num_hq_frag <- args[3] %>% as.numeric()
flag_filter_bp <- args[4] %>% as.character()
cutoff_prop_bp <- args[5] %>% as.numeric()
cutoff_sd_indiv <- args[6] %>% as.numeric()
cutoff_mad_indiv <- args[7] %>% as.numeric()
cutoff_sd_both <- args[8] %>% as.numeric()
cutoff_mad_both <- args[9] %>% as.numeric()
cutoff_sd_extreme <- args[10] %>% as.numeric()
cutoff_mad_extreme <- args[11] %>% as.numeric()

files <- list.files(pattern = "res_table*")

df <- NULL
for (mfile in files) {
  df <- read.table(mfile, header = TRUE, sep = "\t", comment.char = "", quote = "") %>%
    rbind(df, .)
}

# Order table and then filter
mchr <- c(
  "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
  "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
)

df$CHROM <- df$CHROM %>% factor(levels = mchr)

df <- with(df, order(CHROM, POS)) %>% df[., ]

df <- df %>% dplyr::select(.data = ., -"IdRow")

# Write raw table
write.table(df, file = "res_raw.tsv", append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

### Filter table based on cutoffs

df_a <- df %>%
  subset(MEDIAN_AS_VARIANT_READS >= cutoff_as) %>%
  droplevels()

write.table(df_a, file = "res_1as.tsv", append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

df_b <- df_a %>%
  subset(PropClippedSupportVariant < cutoff_clip) %>%
  droplevels()

write.table(df_b, file = "res_1as_2propclip.tsv", append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

df_selected <- df_b %>%
  subset(NUM_FRAGMENTS_MQ40_BQ30_FR > cutoff_num_hq_frag) %>%
  droplevels()

write.table(df_selected, file = "res_1as_2propclip_3numhq.tsv", append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


if (flag_filter_bp == "true") {
  ### Identify low filtering in each strand
  df_a <- df_selected %>%
    subset(NUM_FRAGMENTS_MQ40_BQ30_F < 2 & NUM_FRAGMENTS_MQ40_BQ30_R > 1) %>%
    subset((PROP_FRAGMENTS_BPSTART_UNDER_MQ40_BQ30_R <= cutoff_prop_bp) | (SD_BPSTART_FRAGMENTS_MQ40_BQ30_R > cutoff_sd_indiv & MAD_BPSTART_RAGMENTS_MQ40_BQ30_R > cutoff_mad_indiv)) %>%
    droplevels()

  df_b <- df_selected %>%
    subset(NUM_FRAGMENTS_MQ40_BQ30_F > 1 & NUM_FRAGMENTS_MQ40_BQ30_R < 2) %>%
    subset((PROP_FRAGMENTS_BPSTART_UNDER_MQ40_BQ30_F <= cutoff_prop_bp) | (SD_BPSTART_FRAGMENTS_MQ40_BQ30_F > cutoff_sd_indiv & MAD_BPSTART_RAGMENTS_MQ40_BQ30_F > cutoff_mad_indiv)) %>%
    droplevels()

  df_c <- df_selected %>%
    subset(NUM_FRAGMENTS_MQ40_BQ30_F > 1 & NUM_FRAGMENTS_MQ40_BQ30_R > 1)

  i_a <- df_c %>%
    subset(((PROP_FRAGMENTS_BPSTART_UNDER_MQ40_BQ30_F <= cutoff_prop_bp) & (SD_BPSTART_FRAGMENTS_MQ40_BQ30_F > cutoff_sd_both & MAD_BPSTART_RAGMENTS_MQ40_BQ30_F > cutoff_mad_both)) | (SD_BPSTART_FRAGMENTS_MQ40_BQ30_R > cutoff_sd_extreme & MAD_BPSTART_RAGMENTS_MQ40_BQ30_R > cutoff_mad_extreme)) %>%
    rownames()

  i_b <- df_c %>%
    subset(((PROP_FRAGMENTS_BPSTART_UNDER_MQ40_BQ30_R <= cutoff_prop_bp) & (SD_BPSTART_FRAGMENTS_MQ40_BQ30_R > cutoff_sd_both & MAD_BPSTART_RAGMENTS_MQ40_BQ30_R > cutoff_mad_both)) | (SD_BPSTART_FRAGMENTS_MQ40_BQ30_F > cutoff_sd_extreme & MAD_BPSTART_RAGMENTS_MQ40_BQ30_F > cutoff_mad_extreme)) %>%
    rownames()

  i_end <- intersect(i_a, i_b)

  df_c <- df_c[which(rownames(df_c) %in% i_end), ]

  df_end <- rbind(df_a, df_b, df_c)

  df_end <- with(df_end, order(CHROM, POS)) %>% df_end[., ]
} else {
  df_end <- df_selected
}

## Write filtered table
write.table(df_end, file = "res_1as_2propclip_3numhq_4posvar.tsv", append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
