# INSTALLING AND LOADING PACKAGES
packages <- c("readxl", "dplyr", "tidyr", "ggplot2", "writexl")
install.packages(packages[!packages %in% installed.packages()[, "Package"]])
lapply(packages, library, character.only = TRUE)

# SET WORKING DIRECTORY
setwd("~/Documents/MSF Data/NR-AET Project/Github Upload/")
corr_input <- read_excel("NR-AET_corr matrix data.xlsx", sheet = 1)


# SELECT WHOLE-BODY (B-H) AND GO PATHWAY (I-R) COLUMNS
whole_body <- corr_input[, 2:8]
go_pathways <- corr_input[, 9:18]

# FORCE NUMERIC IN CASE EXCEL TYPES ARE MIXED
whole_body <- whole_body %>% mutate(across(everything(), as.numeric))
go_pathways <- go_pathways %>% mutate(across(everything(), as.numeric))

# PEARSON CORRELATION MATRIX: WHOLE-BODY VS GO
cor_mat <- cor(whole_body, go_pathways, use = "pairwise.complete.obs", method = "pearson")

# FORMAT LABELS TO MATCH FIGURE STYLE
row_labels <- c(
  "Body Weight",
  "Fasting\n[Blood Glucose]",
  "Light Cycle RER",
  "Liver\n% Lipid Area",
  "Liver TAG\nAbundance",
  "Liver DAG\nAbundance",
  "Liver Ceramide\nAbundance"
)

col_labels <- paste0(colnames(go_pathways), " (", seq_len(ncol(go_pathways)), ")")

plot_df <- as.data.frame(as.table(cor_mat)) %>%
  rename(WholeBody = Var1, GO = Var2, r = Freq) %>%
  mutate(
    WholeBody = factor(WholeBody, levels = rownames(cor_mat), labels = row_labels),
    GO = factor(GO, levels = colnames(cor_mat), labels = col_labels),
    abs_r = abs(r)
  )

# DRAW BUBBLE CORRELATION MATRIX
corr_plot <- ggplot(plot_df, aes(x = GO, y = WholeBody)) +
  geom_tile(fill = "white", color = "#b3b3b3", linewidth = 0.55) +
  geom_point(aes(size = abs_r, fill = r), shape = 21, color = "black", stroke = 0) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, 0, 1),
    name = NULL
  ) +
  scale_size_continuous(
    range = c(1, 18),
    limits = c(0, 1),
    breaks = c(0.3, 0.5, 0.7),
    labels = c("r = 0.3", "r = 0.5", "r = 0.7"),
    name = NULL
  ) +
  scale_y_discrete(limits = rev(levels(plot_df$WholeBody))) +
  scale_x_discrete(position = "top") +
  guides(
    fill = guide_colorbar(order = 1, barheight = unit(5.8, "cm")),
    size = guide_legend(order = 2, override.aes = list(fill = NA, color = "black", shape = 21))
  ) +
  labs(x = NULL, y = NULL, tag = "B") +
  theme_minimal(base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, angle = 25, hjust = 0, vjust = 0.1, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 11, face = "bold"),
    legend.key = element_rect(fill = NA, color = NA),
    plot.tag = element_text(size = 48, face = "bold"),
    plot.tag.position = c(0.01, 0.98),
    plot.margin = margin(20, 10, 10, 10)
  )

# SAVE OUTPUTS
print(corr_plot)
ggsave("NR-AET_correlation_matrix_plot.png", plot = corr_plot, width = 13, height = 8, units = "in", dpi = 300)
cor_out <- data.frame(
  WholeBody_Metric = rownames(cor_mat),
  as.data.frame(cor_mat),
  row.names = NULL,
  check.names = FALSE
)