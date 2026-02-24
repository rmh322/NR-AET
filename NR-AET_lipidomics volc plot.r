# INSTALLING AND LOADING PACKAGES
packages <- c("ggplot2", "dplyr", "openxlsx", "readxl", "ggrepel", "agricolae")
install.packages(packages[!packages %in% installed.packages()[, "Package"]])
invisible(lapply(packages, library, character.only = TRUE))

# SETTING WORKING DIRECTORY
setwd("~/Documents/MSF Data/NR-AET Project/Github Upload/")
data_clean <- read_excel("NR-AET_lipidomics volc plot raw data.xlsx", sheet = 1)

# STATS REQUIRED FOR VOLCANO PLOTS
p.values <- numeric(nrow(data_clean))
fc_aet_vec <- numeric(nrow(data_clean)) # HFD+AET vs HFD fold change
fc_nr_vec <- numeric(nrow(data_clean))  # HFD+AET+NR vs HFD fold change
posthoc_results <- vector("list", nrow(data_clean))

for (i in seq_len(nrow(data_clean))) {
  current_data <- data.frame(
    Expression = as.numeric(as.matrix(data_clean[i, -1])),
    Group = factor(c(rep("HFD", 9), rep("HFD+AET", 10), rep("HFD+AET+NR", 10)))
  )

  aov_result <- aov(Expression ~ Group, data = current_data)
  anova_summary <- summary(aov_result)[[1]]
  p.values[i] <- anova_summary$`Pr(>F)`[1]

  # Kept for exported workbook traceability with previous workflow.
  lsd_result <- LSD.test(aov_result, "Group", alpha = 0.5, group = TRUE)
  posthoc_results[[i]] <- lsd_result$groups

  group_means <- tapply(current_data$Expression, current_data$Group, mean)
  fc_aet_vec[i] <- group_means["HFD+AET"] / group_means["HFD"]
  fc_nr_vec[i] <- group_means["HFD+AET+NR"] / group_means["HFD"]
}

p.fdr <- p.adjust(p.values, method = "fdr")
lsd_groups <- sapply(posthoc_results, function(x) paste(rownames(x), x$groups, collapse = "; "))

# EXPORT: volcano input values
comprehensive_results <- data.frame(
  Lipid = data_clean$Lipid,
  p_value = p.values,
  p_adj = p.fdr,
  fc_AET_HFD = fc_aet_vec,
  log2_fc_AET_HFD = log2(fc_aet_vec),
  fc_NR_HFD = fc_nr_vec,
  log2_fc_NR_HFD = log2(fc_nr_vec),
  lsd_groups = lsd_groups,
  is_significant = p.values < 0.05,
  is_DAG = startsWith(data_clean$Lipid, "d_")
)

fc_threshold <- 0.5
p_threshold <- 0.05

comprehensive_results$significant_in_volcano_AET <-
  abs(log2(comprehensive_results$fc_AET_HFD)) >= fc_threshold &
  comprehensive_results$p_value < p_threshold

comprehensive_results$significant_in_volcano_NR <-
  abs(log2(comprehensive_results$fc_NR_HFD)) >= fc_threshold &
  comprehensive_results$p_value < p_threshold

write.xlsx(
  comprehensive_results,
  "NR-AET_lipidomics volc plot values.xlsx",
  rowNames = FALSE,
  colNames = TRUE
)

# VOLCANO PLOT DATA
volcano_data <- data.frame(
  Lipid = rep(data_clean$Lipid, 2),
  log2FC = c(log2(fc_aet_vec), log2(fc_nr_vec)),
  logP = rep(-log10(p.values), 2),
  Comparison = rep(c("HFD+AET", "HFD+AET+NR"), each = nrow(data_clean))
)

volcano_data$is_dag <- startsWith(volcano_data$Lipid, "d_")
volcano_data$significance <- case_when(
  abs(volcano_data$log2FC) >= fc_threshold & volcano_data$logP > -log10(p_threshold) ~ "Significant",
  TRUE ~ "Not Significant"
)

volcano_data$dag_category <- case_when(
  volcano_data$is_dag & volcano_data$significance == "Significant" ~ "Significant DAG",
  volcano_data$is_dag ~ "Non-significant DAG",
  volcano_data$significance == "Significant" ~ "Significant Other",
  TRUE ~ "Not Significant"
)

# VOLCANO PLOT (highlighting significant DAG)
volcano_plot <- ggplot(volcano_data, aes(x = log2FC, y = logP)) +
  geom_point(
    data = subset(volcano_data, dag_category == "Not Significant"),
    size = 2, alpha = 0.3, color = "grey90"
  ) +
  geom_point(
    data = subset(volcano_data, dag_category == "Significant Other"),
    aes(color = Comparison), size = 3, alpha = 0.5
  ) +
  geom_point(
    data = subset(volcano_data, dag_category == "Non-significant DAG"),
    size = 3, alpha = 0.6, color = "grey90"
  ) +
  geom_point(
    data = subset(volcano_data, dag_category == "Significant DAG"),
    aes(color = Comparison), size = 3, alpha = 1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey70") +
  scale_color_manual(
    values = c("HFD+AET" = "grey60", "HFD+AET+NR" = "#8698ff"),
    labels = c("HFD+AET vs HFD", "HFD+AET+NR vs HFD")
  ) +
  scale_x_continuous(breaks = seq(-3, 1, 1)) +
  scale_y_continuous(breaks = seq(0, 8, 2)) +
  labs(
    title = "Liver: Lipid Species Volcano Plot",
    x = expression("log"[2] * "FC"),
    y = expression("-log"[10] * "p-value"),
    color = "Comparison"
  ) +
  geom_text_repel(
    data = subset(volcano_data, dag_category == "Significant DAG"),
    aes(label = Lipid, color = Comparison),
    size = 3,
    fontface = "bold",
    max.overlaps = 30,
    box.padding = 0.5,
    min.segment.length = 0,
    seed = 42
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(1.2, 1.2, 1.2, 1.2, "cm"),
    aspect.ratio = 0.8
  )

print(volcano_plot)

ggsave(
  "NR-AET_lipidomics volc plot.pdf",
  volcano_plot,
  width = 7,
  height = 7,
  dpi = 300,
  units = "in"
)

# VOLCANO PLOT (all labels)
all_labeled_volcano <- ggplot(volcano_data, aes(x = log2FC, y = logP)) +
  geom_point(aes(color = Comparison, alpha = significance), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey70") +
  scale_color_manual(
    values = c("HFD+AET" = "grey60", "HFD+AET+NR" = "#8698ff"),
    labels = c("HFD+AET vs HFD", "HFD+AET+NR vs HFD")
  ) +
  scale_alpha_manual(values = c("Significant" = 1, "Not Significant" = 0.3)) +
  scale_x_continuous(breaks = seq(-3, 1, 1)) +
  scale_y_continuous(breaks = seq(0, 8, 2)) +
  labs(
    title = "Liver: All Labeled Lipid Species",
    x = expression("log"[2] * "FC"),
    y = expression("-log"[10] * "p-value"),
    color = "Comparison"
  ) +
  geom_text_repel(
    aes(label = Lipid, color = Comparison),
    size = 2.5,
    max.overlaps = 100,
    box.padding = 0.3,
    min.segment.length = 0,
    seed = 42
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(1.2, 1.2, 1.2, 1.2, "cm"),
    aspect.ratio = 0.8
  ) +
  guides(alpha = "none")

print(all_labeled_volcano)

ggsave(
  "NR-AET_lipidomics volc plot (all labels).pdf",
  all_labeled_volcano,
  width = 9,
  height = 8,
  dpi = 300,
  units = "in"
)
