# Correlation matrix bubble plot: FBG (column B) vs DAG species (columns C:O)
# Input: NR-AET_lipidomics FBG v DAG.xlsx (sheet "DAG Only")

packages <- c("readxl", "dplyr", "tidyr", "ggplot2", "grid", "writexl")
missing_packages <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
invisible(lapply(packages, library, character.only = TRUE))

input_file <- "NR-AET_lipidomics FBG v DAG.xlsx"
input_sheet <- "DAG Only"
output_plot <- "NR-AET_lipidomics FBG v DAG_corr matrix.png"
output_values_xlsx <- "NR-AET_lipidomics FBG v DAG_corr matrix_values.xlsx"

raw_df <- read_excel(input_file, sheet = input_sheet)

fbg_col <- names(raw_df)[2]
dag_cols <- names(raw_df)[3:15]

analysis_df <- raw_df %>%
  mutate(
    GroupCode = case_when(
      grepl("^HFD", Group) ~ "HFD",
      grepl("^EX", Group) ~ "EX",
      grepl("^NR", Group) ~ "NR",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(GroupCode)) %>%
  mutate(across(all_of(c(fbg_col, dag_cols)), as.numeric))

# Group order chosen to reproduce the reference panel layout.
group_levels <- c("HFD", "EX", "NR")
symbol_labels <- c("HFD" = "-\n-", "EX" = "+\n-", "NR" = "+\n+")

corr_df <- expand.grid(GroupCode = group_levels, DAG = dag_cols, stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(
    cor_test = list(
      cor.test(
        x = analysis_df[[fbg_col]][analysis_df$GroupCode == GroupCode],
        y = analysis_df[[DAG]][analysis_df$GroupCode == GroupCode],
        method = "pearson"
      )
    ),
    r = unname(cor_test$estimate),
    p_value = cor_test$p.value,
    abs_r = abs(r)
  ) %>%
  mutate(p_value_fdr = p.adjust(p_value, method = "fdr")) %>%
  select(-cor_test) %>%
  ungroup() %>%
  mutate(
    GroupCode = factor(GroupCode, levels = group_levels),
    DAG = factor(DAG, levels = rev(dag_cols))
  )

write_xlsx(
  corr_df %>% select(GroupCode, DAG, r, p_value, p_value_fdr),
  output_values_xlsx
)

n_dag <- length(dag_cols)

p <- ggplot(corr_df, aes(x = GroupCode, y = DAG)) +
  geom_tile(fill = "white", color = "grey70", linewidth = 0.7) +
  geom_point(aes(size = abs_r, fill = r), shape = 21, color = "black", stroke = 0) +
  scale_fill_gradient2(
    low = "#1f27ff",
    mid = "white",
    high = "#ff2a2a",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, 0, 1),
    name = NULL
  ) +
  scale_size_continuous(
    range = c(1.5, 11),
    limits = c(0, 1),
    breaks = c(0.3, 0.5, 0.7),
    labels = c("r = 0.3", "r = 0.5", "r = 0.7"),
    name = NULL
  ) +
  scale_x_discrete(position = "top", labels = symbol_labels) +
  guides(
    fill = guide_colorbar(order = 1, barheight = unit(3, "cm"), barwidth = unit(0.7, "cm")),
    size = guide_legend(
      order = 2,
      override.aes = list(fill = "white", color = "black", shape = 21, stroke = 0.8)
    )
  ) +
  labs(
    title = expression(italic("Fasting Blood Glucose (mmol" %.% L^{-1} * ")")),
    x = NULL,
    y = NULL,
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    plot.background = element_rect(fill = "grey68", color = NA),
    axis.text.x = element_text(size = 15, face = "bold", color = "black", lineheight = 1.1),
    axis.text.y = element_text(size = 13, face = "italic", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 16)),
    legend.text = element_text(size = 13, face = "bold"),
    legend.key = element_rect(fill = "grey68", color = NA),
    legend.background = element_rect(fill = "grey68", color = NA),
    legend.box.background = element_rect(fill = "grey68", color = NA),
    plot.title = element_text(size = 17, face = "bold.italic", hjust = 0.52, margin = margin(b = 8)),
    plot.tag = element_text(size = 34, face = "bold"),
    plot.tag.position = c(0.03, 0.98),
    plot.margin = margin(18, 22, 18, 22)
  )

print(p)
ggsave(output_plot, p, width = 4.7, height = 9.4, units = "in", dpi = 350, bg = "grey68")
