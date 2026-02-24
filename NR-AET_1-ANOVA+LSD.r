# INSTALLING AND LOADING PACKAGES 
packages <- c("vegan", "ggplot2", "dplyr", "colorRamps", "pheatmap", "openxlsx", "readxl", "tidyr",
              "writexl", "tibble", "ggrepel", "agricolae")
install.packages(packages[!packages %in% installed.packages()[,"Package"]])
lapply(packages, library, character.only = TRUE)
library(agricolae)  # Explicitly load agricolae

# SET WORKING DIRECTORY 
setwd("~/Documents/MSF Data/NR-AET Project/Github Upload/")
data_clean <- read_excel("NR-AET_raw data.xlsx", sheet = 2) # Want to use the clean version here

# RESHAPE DATA FOR ANALYSIS
long_data <- data_clean %>%
  pivot_longer(cols = -Protein, names_to = "Sample", values_to = "Expression") %>%
  mutate(
    Group = case_when(
      grepl("^H", Sample) ~ "HFD",
      grepl("^E", Sample) ~ "HFD+AET",
      grepl("^N", Sample) ~ "HFD+AET+NR"
    )
  )
table(long_data$Group, long_data$Sample)

# PERFORM ONE-WAY ANOVA AND FISHER'S LSD
analyze_protein <- function(df) {
  # One-way ANOVA
  model <- aov(Expression ~ Group, data = df)
  anova_summary <- summary(model)[[1]]
  
  # Fisher's LSD test
  lsd_result <- LSD.test(model, "Group", p.adj = "fdr")
  
  # Extract pairwise comparisons
  comparisons <- as.data.frame(lsd_result$groups)
  
  # Calculate fold changes between all groups
  group_means <- tapply(df$Expression, df$Group, mean)
  fc_combinations <- combn(names(group_means), 2, simplify = FALSE)
  
  fc_results <- sapply(fc_combinations, function(pair) {
    group_means[pair[2]] / group_means[pair[1]]
  })
  
  # Create names for fold changes
  fc_names <- sapply(fc_combinations, function(pair) {
    paste0("fc_", pair[2], "_vs_", pair[1])
  })
  
  # Combine results
  result <- data.frame(
    Protein = df$Protein[1],
    p_ANOVA = anova_summary["Group", "Pr(>F)"]
  )
  
  # Add fold changes to results
  for(i in seq_along(fc_names)) {
    result[[fc_names[i]]] <- fc_results[i]
  }
  
  # Add LSD p-values
  pairwise_p <- lsd_result$groups
  result$lsd_groups <- paste(rownames(pairwise_p), pairwise_p$groups, sep = ":", collapse = "; ")
  
  return(result)
}

# APPLY THE ANALYSIS FOR EACH PROTEIN 
results <- long_data %>%
  group_by(Protein) %>%
  do(analyze_protein(.))

# ADJUST ANOVA p-values FOR MULTIPLE COMPARISONS
results <- results %>%
  mutate(
    p_ANOVA_adj = p.adjust(p_ANOVA, method = "fdr")
  )

# CREATE A DATA FRAME WITH SIGNIFICANT PROTEINS ONLY
significant_proteins <- results %>%
  filter(p_ANOVA_adj < 0.05)

# GENERATE DATA FOR VISUALIZATION
visualization_data <- results %>%
  mutate(
    across(starts_with("fc_"), log2),
    pvolc = -log10(p_ANOVA_adj)
  )

# SAVE RESULTS
write_xlsx(visualization_data, "NR-AET_full analysis.xlsx")
write_xlsx(significant_proteins, "NR-AET_sig pr.xlsx")
