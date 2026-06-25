# NR-AET
This repository contains supplementary materials, raw lipidomics and proteomics data and R scripts for data processing and analysis used in the study titled: **"Nitrate and resveratrol supplementation selectively enhances hepatic adaptations to aerobic exercise in high-fat fed male mice"**.

## Files
**SUPPLEMENTARY DATA**
Supplementary Table 1. Food intake over the intervention period.

Supplementary Figure 1. Comparison of glucose tolerance responses in the current study with reference data from Handy et al. 2023. Fasting blood glucose (A), ipGTT trace (B), and incremental area under the curve (C). Red dashed lines represent reference data from HFD+NR sedentary animals reported previously (ref. #22). AUC, area under the curve; +NR, co-supplementation with dietary nitrate+resveratrol. Data expressed as mean ± SD. * p < 0.05 versus sedentary HFD controls. 

Supplementary Figure 2. Analysis of individual lipid species from lipidomics. iWAT TAG (A), DAG (B) and ceramide (C) content analyzed as a percent change from sedentary HFD controls. eWAT TAG (D), DAG (E) and ceramide (F) content. Skeletal muscle TAG (G), DAG (H) and ceramide (I) content. Liver TAG (J), DAG (K) and ceramide (L) content. DAG, diacyl glycerol; eWAT, epididymal (visceral) white adipose tissue; HFD, high fat diet; iWAT, inguinal (subcutaneous) white adipose tissue; +NR, co-supplementation with dietary nitrate+resveratrol; TAG, triacyl glycerol. Data expressed as mean ± SD. * p < 0.05 versus sedentary HFD controls. ∆ p < 0.05 verses HFD+AET.

**LIPIDOMICS**
- 'NR-AET_lipidomics raw data.xlsx' - Raw lipidomics data
- 'NR-AET_lipidomics FBG v DAG.xlsx' - Raw DAG lipidome and fasting blood glucose data for correlation analysis
- 'NR-AET_corr matrix (FBG v DAG) plot.r' - Script for correlation analysis between DAG and fasting blood glucose levels
- 'NR-AET_lipidomics volc plot.r' - Script for volcano plot

**PROTEOMICS**
- 'NR-AET_raw data.xlsx' - Raw proteomics data
- 'NR-AET_corr matrix data.xlsx' - Raw proteome + whole-body characterization data for correlation analysis
- 'NR-AET_1-ANOVA+LSD.r' - Script for statistical analysis (one-way ANOVA)
- 'NR-AET_PCA+biplot.r' - Script for PCA and biplot
- 'NR_AET_corr matrix plot.r' - Script for correlation anlaysis between whole body characterization parameters and GO pathway protein expression
