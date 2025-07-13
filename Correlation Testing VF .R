#########################################
# Setup: load necessary libraries
#########################################

library(tidyverse)  # Data manipulation and plotting (includes dplyr and ggplot2)
library(janitor)    # Clean column names and text data
library(ggplot2)    # Plotting

#########################################
# Read raw data files
#########################################

# Valley Fever data: skip first 3 header rows and set column names manually
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

# Social Vulnerability Index (SVI) data
svi <- read_csv("California_county.csv")

#########################################
# Clean Valley Fever data
#########################################

vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),                  # Make county names uppercase for consistency
    county = str_remove(county, " COUNTY"),         # Remove " COUNTY" suffix if present
    county = str_trim(county),                      # Remove extra whitespace
    cases = as.numeric(cases_raw),                  # Convert raw cases to numeric
    county = case_when(                             # Reassign certain cities to main county
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%  # Remove summary and unknown rows
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  filter(total_cases > 44)  # Only keep counties with >44 cases (to focus on higher burden areas)

#########################################
# Clean SVI data (standardize county names)
#########################################

svi_clean <- svi %>%
  janitor::clean_names() %>%
  mutate(
    county = str_to_upper(county),                  # Uppercase for consistency
    county = str_remove(county, " COUNTY"),         # Remove " COUNTY" suffix
    county = str_trim(county)                       # Remove extra spaces
  ) %>%
  filter(!is.na(rpl_themes), !is.na(e_totpop)) %>% # Keep rows with non-missing SVI and population
  dplyr::select(county, rpl_themes, e_totpop)      # Select only relevant columns

#########################################
# Check for counties not matching before merging
#########################################

print("Counties in VF data without match in SVI data:")
print(anti_join(vf_clean, svi_clean, by = "county"))

#########################################
# Merge VF and SVI data, then calculate VF rate
#########################################

vf_svi <- left_join(vf_clean, svi_clean, by = "county") %>%
  mutate(vf_rate = (total_cases / e_totpop) * 100000) %>% # Calculate rate per 100,000 people
  drop_na()

#########################################
# Check merged dataset summary (quick sanity check)
#########################################

print(dim(vf_svi))     # Print dimensions of final merged dataset
print(summary(vf_svi)) # Print summary stats

#########################################
# Correlation test between SVI score and VF rate
#########################################

cor_test <- cor.test(vf_svi$rpl_themes, vf_svi$vf_rate)

# Save Pearson correlation coefficient and p-value
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)

# Print correlation results to console
cat("Pearson correlation coefficient (r):", r_value, "\n")
cat("p-value:", p_value, "\n")

# Build annotation text for graph
annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

#########################################
# Plot: VF rate vs. SVI score
#########################################

ggplot(vf_svi, aes(x = rpl_themes, y = vf_rate)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate("text", 
           x = min(vf_svi$rpl_themes, na.rm = TRUE) + 0.02,
           y = max(vf_svi$vf_rate, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4, color = "black") +
  labs(
    title = "Valley Fever Rate vs. County-Level SVI Overall Score",
    x = "County-Level SVI Overall Score (rpl_themes)",
    y = "VF Rate per 100,000 (2001–2023)"
  ) +
  theme_minimal()

#########################################
# Linear regression: VF rate ~ SVI score
#########################################

model_lm <- lm(vf_rate ~ rpl_themes, data = vf_svi)

# Print detailed regression summary (coefficients, R², p-values, etc.)
summary(model_lm)

#########################################
# Additional context for teammate
#########################################

# What is a good correlation?
# Generally, |r| > 0.5 = moderate to strong, but context matters.
# Our r value suggests a moderate correlation but may not be strong enough to rely on alone.

# Why we use linear regression?
# To quantify how much VF rate changes as SVI score changes, not just if they are related.

# Why we also report p-value?
# To confirm that the relationship is statistically significant and not random.

# Plot residuals
par(mfrow = c(2, 2))  # Diagnostic plot panel
plot(model_lm)
par(mfrow = c(1, 1))  # Reset


# Null model (just mean)
model_null <- lm(vf_rate ~ 1, data = vf_svi)

# Compare with AIC
AIC(model_lm, model_null)



#########################################
# Setup: load necessary libraries
#########################################

library(tidyverse)  # Includes dplyr and ggplot2, for data wrangling and plotting
library(janitor)    # For cleaning column names
library(corrplot)   # For creating correlation heatmaps
library(car)        # For calculating Variance Inflation Factor (VIF)

#########################################
# Read and clean SVI data
#########################################

svi <- read_csv("California_county.csv") %>%
  janitor::clean_names() %>%                     # Convert all column names to snake_case
  mutate(
    county = str_to_upper(county),              # Standardize county names to uppercase
    county = str_remove(county, " COUNTY"),     # Remove " COUNTY" suffix
    county = str_trim(county)                   # Remove extra whitespace
  )

#########################################
# Select 16 SVI component variables
#########################################

svi_vars <- svi %>%
  dplyr::select(
    ep_pov150, ep_unemp, ep_hburd, ep_nohsdp, ep_uninsur,
    ep_age65, ep_age17, ep_disabl, ep_sngpnt, ep_limeng,
    ep_minrty, ep_munit, ep_mobile, ep_crowd, ep_noveh, ep_groupq
  ) %>%
  drop_na()  # Remove rows with missing values to ensure clean correlation calculations

#########################################
# Generate and plot correlation matrix
#########################################

# Calculate pairwise correlations using complete observations
corr_matrix <- cor(svi_vars, use = "complete.obs")

# Clear any existing plots
dev.off()

# Plot heatmap of correlations
corrplot(corr_matrix, method = "color", type = "upper",
         tl.col = "black",             # Text label color
         tl.cex = 0.7,                 # Text label size
         addCoef.col = "black",        # Add correlation coefficient text
         number.cex = 0.6,             # Size of coefficient text
         title = "Correlation Heat Map of 16 SVI Variables",
         mar = c(0, 0, 2, 0))          # Margins for the title

#########################################
# Calculate Variance Inflation Factors (VIF)
#########################################

# Create dummy response variable since we only want to examine multicollinearity among predictors
dummy_response <- rnorm(nrow(svi_vars))
model_vif <- lm(dummy_response ~ ., data = svi_vars)

# Calculate VIF values
vif_values <- vif(model_vif)

#########################################
# Create friendly (readable) variable names for plot
#########################################

friendly_names <- c(
  ep_pov150 = "Poverty",
  ep_unemp = "Unemployment",
  ep_hburd = "High Housing Burden",
  ep_nohsdp = "No HS Diploma",
  ep_uninsur = "Uninsured",
  ep_age65 = "Aged 65+",
  ep_age17 = "Under 17",
  ep_disabl = "Disability",
  ep_sngpnt = "Single Parent",
  ep_limeng = "Limited English",
  ep_minrty = "Minority Status",
  ep_munit = "Multi-Unit Housing",
  ep_mobile = "Mobile Homes",
  ep_crowd = "Crowding",
  ep_noveh = "No Vehicle",
  ep_groupq = "Group Quarters"
)

#########################################
# Create VIF data frame for plotting
#########################################

vif_df <- data.frame(Variable = names(vif_values), VIF = vif_values) %>%
  mutate(Friendly = friendly_names[Variable])  # Map to friendly names

#########################################
# Plot VIFs
#########################################

ggplot(vif_df, aes(x = reorder(Friendly, VIF), y = VIF)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip axes for easier reading
  geom_hline(yintercept = 5, color = "orange", linetype = "dashed", linewidth = 1.5) +
  annotate("text", x = 1, y = 5.2, label = "VIF = 5", color = "orange", hjust = 0, size = 5, fontface = "bold") +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed", linewidth = 1.5) +
  annotate("text", x = 1, y = 10.2, label = "VIF = 10", color = "red", hjust = 0, size = 5, fontface = "bold") +
  labs(
    title = "Variance Inflation Factors for SVI Variables",
    x = "SVI Component",
    y = "VIF"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

#########################################
# Explanation for teammate
#########################################

# In this code, we:
# 1️⃣ Clean and prepare SVI data.
# 2️⃣ Calculate correlations among all 16 SVI variables and plot a heatmap to check for strong linear relationships.
# 3️⃣ Calculate Variance Inflation Factors (VIF) to detect multicollinearity among predictors — high VIF means a variable is highly redundant.
# 4️⃣ Visualize VIFs using a bar plot with threshold lines at 5 and 10 for easy interpretation.








#########################################
# Setup: Load libraries
#########################################

install.packages("leaps")  # Only needed once — comment out after first install
library(leaps)             # For best subset regression
library(car)               # For checking multicollinearity with VIF
library(tidyverse)         # Includes dplyr, ggplot2, etc.
library(janitor)           # For cleaning column names

#########################################
# Read and clean SVI data
#########################################

svi_all <- read_csv("data_raw/california_svi_2020.csv") %>%
  clean_names() %>%                           # Make all column names lowercase, underscores
  mutate(county = str_to_upper(county))       # Standardize county name case

#########################################
# Summarize SVI variables by county
#########################################
# Here we take mean of each vulnerability indicator per county

svi_summary <- svi_all %>%
  group_by(county) %>%
  summarise(
    pov = mean(ep_pov150, na.rm = TRUE),
    unemp = mean(ep_unemp, na.rm = TRUE),
    hburd = mean(ep_hburd, na.rm = TRUE),
    nohsdp = mean(ep_nohsdp, na.rm = TRUE),
    uninsur = mean(ep_uninsur, na.rm = TRUE),
    age65 = mean(ep_age65, na.rm = TRUE),
    age17 = mean(ep_age17, na.rm = TRUE),
    disabl = mean(ep_disabl, na.rm = TRUE),
    sngpnt = mean(ep_sngpnt, na.rm = TRUE),
    limeng = mean(ep_limeng, na.rm = TRUE),
    minrty = mean(ep_minrty, na.rm = TRUE),
    munit = mean(ep_munit, na.rm = TRUE),
    mobile = mean(ep_mobile, na.rm = TRUE),
    crowd = mean(ep_crowd, na.rm = TRUE),
    noveh = mean(ep_noveh, na.rm = TRUE),
    groupq = mean(ep_groupq, na.rm = TRUE),
    .groups = "drop"
  )

#########################################
# Read and clean Valley Fever data
#########################################

vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),                   # Standardize case
    county = str_remove(county, " COUNTY"),          # Remove " COUNTY" text
    county = str_trim(county),                       # Trim spaces
    cases = as.numeric(cases_raw),                   # Convert cases column to numeric
    county = case_when(                              # Combine smaller jurisdictions into LA
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop")

#########################################
# Merge SVI and Valley Fever data
#########################################

vf_svi <- left_join(vf_clean, svi_summary, by = "county") %>%
  drop_na()  # Remove counties with missing values in any variable

#########################################
# Run best subset regression
#########################################

# Create list of predictor variable names
predictor_vars <- names(svi_summary)[-1]  # Remove "county" column

# Create formula text
formula_text <- paste("total_cases ~", paste(predictor_vars, collapse = " + "))

# Run best subset regression (tries all variable combinations up to 16)
best_subset <- regsubsets(as.formula(formula_text), data = vf_svi, nvmax = 16, really.big = TRUE)
summary_best <- summary(best_subset)

#########################################
# Evaluate subsets and check VIF
#########################################

selected_matrix <- summary_best$outmat
successful_combos <- list()

for (size in 1:16) {
  # Find combinations using 'size' number of variables
  combo_rows <- rownames(selected_matrix)[apply(selected_matrix, 1, function(x) sum(x == "*")) == size]
  
  for (row in combo_rows) {
    vars_in_combo <- names(which(selected_matrix[row, ] == "*"))
    formula_combo <- paste("total_cases ~", paste(vars_in_combo, collapse = " + "))
    model_combo <- lm(as.formula(formula_combo), data = vf_svi)
    
    if (length(vars_in_combo) == 1) {
      # Only one predictor → no multicollinearity to check
      adj_r2 <- summary(model_combo)$adj.r.squared
      successful_combos[[length(successful_combos) + 1]] <- list(
        vars = vars_in_combo,
        adj_r2 = adj_r2,
        vif = NA
      )
    } else {
      vifs <- vif(model_combo)
      if (all(vifs < 5)) {  # Only keep combinations with all VIF < 5
        adj_r2 <- summary(model_combo)$adj.r.squared
        successful_combos[[length(successful_combos) + 1]] <- list(
          vars = vars_in_combo,
          adj_r2 = adj_r2,
          vif = vifs
        )
      }
    }
  }
}

#########################################
# Print successful combinations
#########################################

if (length(successful_combos) == 0) {
  cat("No combinations found with all VIF < 5.\n")
} else {
  for (i in seq_along(successful_combos)) {
    cat("\n--- Successful combination", i, "---\n")
    print(successful_combos[[i]]$vars)
    cat("Adjusted R²:", round(successful_combos[[i]]$adj_r2, 3), "\n")
    print(successful_combos[[i]]$vif)
  }
}

#########################################
# Explanation for your teammate
#########################################

# 1️⃣ This code first cleans and summarizes SVI data by county.
# 2️⃣ It then cleans Valley Fever case data and merges with SVI data.
# 3️⃣ We use best subset regression to try every possible combination of SVI predictors (up to 16).
# 4️⃣ For each combination, it checks if all VIF values are below 5 to avoid multicollinearity.
# 5️⃣ Successful models (with VIF < 5) and their adjusted R² values are printed to help select the most stable models.






#########################################
# Load libraries
#########################################
library(tidyverse)
library(janitor)
library(car)

#########################################
# --- Read and clean SVI data ---
#########################################
svi_all <- read_csv("data_raw/california_svi_2020.csv") %>%
  clean_names()

# Summarize SVI variables by county: take mean of selected vulnerability metrics
svi_summary <- svi_all %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    pov = mean(ep_pov150, na.rm = TRUE),
    disabl = mean(ep_disabl, na.rm = TRUE),
    crowd = mean(ep_crowd, na.rm = TRUE),
    .groups = "drop"
  )

#########################################
# --- Read and clean Valley Fever data ---
#########################################
#########################################
# --- Read and clean Valley Fever data ---
#########################################
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),
    county = str_remove(county, " COUNTY"),
    county = str_trim(county),
    cases = as.numeric(cases_raw),
    county = case_when(
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  filter(total_cases > 44)  # ✅ Keep only counties with more than 44 cases

#########################################
# --- Merge VF data and SVI data ---
#########################################
vf_svi <- left_join(vf_clean, svi_summary, by = "county") %>%
  drop_na()

#########################################
# --- Build and summarize multi-variable regression model ---
#########################################
final_formula <- "total_cases ~ pov + disabl + crowd"
model_final <- lm(as.formula(final_formula), data = vf_svi)

# Print summary: coefficients, significance, R², etc.
summary(model_final)

# Check VIFs
vif_values <- vif(model_final)
cat("VIF for each variable:\n")
print(vif_values)

#########################################
# --- Add predicted values and plot ---
#########################################
vf_svi <- vf_svi %>%
  mutate(predicted_cases = predict(model_final, newdata = vf_svi))

# Calculate correlation between predicted and actual
cor_test <- cor.test(vf_svi$predicted_cases, vf_svi$total_cases)
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)

# Explicitly print correlation score
cat("Correlation between predicted and actual cases (r):", r_value, "\n")
cat("p-value:", p_value, "\n")

# Create annotation text
annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

# Plot predicted vs actual
ggplot(vf_svi, aes(x = predicted_cases, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate("text",
           x = min(vf_svi$predicted_cases, na.rm = TRUE) + 0.02,
           y = max(vf_svi$total_cases, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4, color = "black") +
  labs(
    title = "Predicted vs. Actual Total Valley Fever Cases",
    x = "Predicted Total Cases (Fitted Values)",
    y = "Observed Total Cases (2001–2023)"
  ) +
  theme_minimal()

#########################################
# --- Residual diagnostic plots ---
#########################################
par(mfrow = c(2, 2))  # Diagnostic panel
plot(model_final)
par(mfrow = c(1, 1))  # Reset

#########################################
# --- Null model and AIC comparison ---
#########################################
model_null <- lm(total_cases ~ 1, data = vf_svi)
aic_values <- AIC(model_final, model_null)
print(aic_values)


#########################################
# --- PM2.5 analysis ---
#########################################

# Read PM2.5 data
pm_raw <- read_csv("HDPulse_data_export.csv", col_names = "X1")

# Remove first 5 rows (metadata)
pm_trim <- pm_raw[-c(1:6), , drop = FALSE]

# Separate columns into county, fips, and pm25 value
pm_clean <- pm_trim %>%
  separate(X1, into = c("county", "fips", "pm25"), sep = ",") %>%
  mutate(
    county = str_remove(county, " County"),
    county = str_to_upper(county),
    county = str_trim(county),
    pm25 = as.numeric(pm25)
  ) %>%
  filter(!is.na(pm25))

# Merge VF data with PM2.5
merged_data <- inner_join(pm_clean, vf_clean, by = "county") %>%
  drop_na()

print(dim(merged_data))  # Show number of matching rows (counties)

if (nrow(merged_data) > 0) {
  # Correlation test between PM2.5 and total VF cases
  cor_test <- cor.test(merged_data$pm25, merged_data$total_cases)
  r_value <- round(cor_test$estimate, 3)
  p_value <- signif(cor_test$p.value, 3)
  
  annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)
  
  # Plot PM2.5 vs. Valley Fever cases
  print(
    ggplot(merged_data, aes(x = pm25, y = total_cases)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "blue") +
      annotate("text", 
               x = min(merged_data$pm25, na.rm = TRUE) + 0.5,
               y = max(merged_data$total_cases, na.rm = TRUE) * 0.9,
               label = annotation_text, hjust = 0, size = 4, color = "black") +
      labs(
        title = "Total Valley Fever Cases vs. PM2.5 by County",
        x = "PM2.5 (µg/m³)",
        y = "Total VF Cases (2001–2023)"
      ) +
      theme_minimal()
  )
  
  # Linear regression summary for PM2.5 alone
  model_pm <- lm(total_cases ~ pm25, data = merged_data)
  print(summary(model_pm))
}

#########################################
# --- Explanation for your teammate ---
#########################################
# - We clean and merge SVI (social vulnerability) data and Valley Fever data.
# - We build a multi-variable regression model with 4 key SVI predictors (poverty, disability, limited English, crowding).
# - We calculate and print VIFs to check multicollinearity.
# - Separately, we analyze PM2.5 vs. Valley Fever: run correlation test, plot relationship, run simple regression.
# - All plots are wrapped in print() to show automatically in scripts and when running in one go.










#########################################
# Load libraries
#########################################
library(tidyverse)
library(janitor)
library(leaps)
library(car)

#########################################
# --- Read and clean SVI data (16 variables) ---
#########################################
svi_all <- read_csv("data_raw/california_svi_2020.csv") %>%
  clean_names()

# Summarize 16 SVI variables by county
svi_summary <- svi_all %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    ep_pov150 = mean(ep_pov150, na.rm = TRUE),
    ep_unemp = mean(ep_unemp, na.rm = TRUE),
    ep_hburd = mean(ep_hburd, na.rm = TRUE),
    ep_nohsdp = mean(ep_nohsdp, na.rm = TRUE),
    ep_uninsur = mean(ep_uninsur, na.rm = TRUE),
    ep_age65 = mean(ep_age65, na.rm = TRUE),
    ep_age17 = mean(ep_age17, na.rm = TRUE),
    ep_disabl = mean(ep_disabl, na.rm = TRUE),
    ep_sngpnt = mean(ep_sngpnt, na.rm = TRUE),
    ep_limeng = mean(ep_limeng, na.rm = TRUE),
    ep_minrty = mean(ep_minrty, na.rm = TRUE),
    ep_munit = mean(ep_munit, na.rm = TRUE),
    ep_mobile = mean(ep_mobile, na.rm = TRUE),
    ep_crowd = mean(ep_crowd, na.rm = TRUE),
    ep_noveh = mean(ep_noveh, na.rm = TRUE),
    ep_groupq = mean(ep_groupq, na.rm = TRUE),
    .groups = "drop"
  )

#########################################
# --- Read and clean PM2.5 data ---
#########################################
pm_raw <- read_csv("HDPulse_data_export.csv", col_names = "X1")

# Remove first 5 rows (headers, CA, US rows)
pm_trim <- pm_raw[-c(1:6), , drop = FALSE]

# Separate columns
pm_clean <- pm_trim %>%
  separate(X1, into = c("county", "fips", "pm25"), sep = ",") %>%
  mutate(
    county = str_remove(county, " County"),
    county = str_to_upper(county),
    county = str_trim(county),
    pm25 = as.numeric(pm25)
  ) %>%
  filter(!is.na(pm25)) %>%
  dplyr::select(county, pm25)   # ✅ Use dplyr::select explicitly here


#########################################
# --- Merge SVI and PM2.5 data ---
#########################################
svi_pm <- left_join(svi_summary, pm_clean, by = "county") %>%
  drop_na()

#########################################
# --- Read and clean Valley Fever data ---
#########################################
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),
    county = str_remove(county, " COUNTY"),
    county = str_trim(county),
    cases = as.numeric(cases_raw),
    county = case_when(
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop")

#########################################
# --- Merge everything together ---
#########################################
vf_full <- left_join(vf_clean, svi_pm, by = "county") %>%
  drop_na()

#########################################
# --- Run best subset regression (all 17 variables) ---
#########################################
predictor_vars <- names(vf_full)[!(names(vf_full) %in% c("county", "total_cases"))]
formula_text <- paste("total_cases ~", paste(predictor_vars, collapse = " + "))

best_subset <- regsubsets(as.formula(formula_text), data = vf_full, nvmax = length(predictor_vars), really.big = TRUE)
summary_best <- summary(best_subset)

#########################################
# --- Check combinations with VIF < 5 ---
#########################################
selected_matrix <- summary_best$outmat
successful_combos <- list()

for (size in 1:length(predictor_vars)) {
  combo_rows <- rownames(selected_matrix)[apply(selected_matrix, 1, function(x) sum(x == "*")) == size]
  
  for (row in combo_rows) {
    vars_in_combo <- names(which(selected_matrix[row, ] == "*"))
    formula_combo <- paste("total_cases ~", paste(vars_in_combo, collapse = " + "))
    model_combo <- lm(as.formula(formula_combo), data = vf_full)
    
    if (length(vars_in_combo) == 1) {
      adj_r2 <- summary(model_combo)$adj.r.squared
      successful_combos[[length(successful_combos) + 1]] <- list(
        vars = vars_in_combo,
        adj_r2 = adj_r2,
        vif = NA
      )
    } else {
      vifs <- vif(model_combo)
      if (all(vifs < 5)) {
        adj_r2 <- summary(model_combo)$adj.r.squared
        successful_combos[[length(successful_combos) + 1]] <- list(
          vars = vars_in_combo,
          adj_r2 = adj_r2,
          vif = vifs
        )
      }
    }
  }
}

#########################################
# --- Print successful combinations ---
#########################################
if (length(successful_combos) == 0) {
  cat("No combinations found with all VIF < 5.\n")
} else {
  for (i in seq_along(successful_combos)) {
    cat("\n--- Successful combination", i, "---\n")
    print(successful_combos[[i]]$vars)
    cat("Adjusted R²:", round(successful_combos[[i]]$adj_r2, 3), "\n")
    print(successful_combos[[i]]$vif)
  }
}

















library(MASS)  # For glm.nb()

# Negative binomial model: total_cases ~ rpl_themes, with offset for log population
model_nb <- glm.nb(total_cases ~ rpl_themes + offset(log(e_totpop)), data = vf_svi)

# Print summary
summary(model_nb)

# AIC
cat("AIC:", AIC(model_nb), "\n")

# Optionally, Pearson residuals diagnostic plot
res <- residuals(model_nb, type = "pearson")
plot(fitted(model_nb), res,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     main = "Residual plot (Negative Binomial)")
abline(h = 0, col = "red")







library(MASS)  # For negative binomial
library(tidyverse)
library(janitor)

#########################################
# --- Read and clean SVI data ---
#########################################
svi_all <- read_csv("data_raw/california_svi_2020.csv") %>%
  clean_names()

# Summarize SVI variables by county
svi_summary <- svi_all %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    pov = mean(ep_pov150, na.rm = TRUE),
    disabl = mean(ep_disabl, na.rm = TRUE),
    crowd = mean(ep_crowd, na.rm = TRUE),
    .groups = "drop"
  )

#########################################
# --- Read and clean Valley Fever data ---
#########################################
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),
    county = str_remove(county, " COUNTY"),
    county = str_trim(county),
    cases = as.numeric(cases_raw),
    county = case_when(
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  filter(total_cases > 44)

#########################################
# --- Merge VF data and SVI summary ---
#########################################
vf_svi <- left_join(vf_clean, svi_summary, by = "county") %>%
  filter(!is.na(pov) & !is.na(disabl) & !is.na(crowd))  # ✅ Fix: explicitly drop rows with missing predictors

#########################################
# --- Negative binomial model ---
#########################################
model_nb_multi <- glm.nb(total_cases ~ pov + disabl + crowd, data = vf_svi)

# Summary
summary(model_nb_multi)

# Null model for AIC comparison
model_null_nb <- glm.nb(total_cases ~ 1, data = vf_svi)
aic_values <- AIC(model_nb_multi, model_null_nb)
print(aic_values)

#########################################
# --- Residual diagnostics ---
#########################################
res_nb <- residuals(model_nb_multi, type = "pearson")
plot(fitted(model_nb_multi), res_nb,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     main = "Residual plot (Negative Binomial)")
abline(h = 0, col = "red")

#########################################
# --- Predicted vs actual plot ---
#########################################
vf_svi <- vf_svi %>%
  mutate(predicted_cases_nb = predict(model_nb_multi, type = "response"))

# Correlation
cor_test <- cor.test(vf_svi$predicted_cases_nb, vf_svi$total_cases)
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)

annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

ggplot(vf_svi, aes(x = predicted_cases_nb, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  annotate("text",
           x = min(vf_svi$predicted_cases_nb, na.rm = TRUE),
           y = max(vf_svi$total_cases, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4, color = "black") +
  labs(
    title = "Predicted vs. Actual Total Valley Fever Cases (Negative Binomial)",
    x = "Predicted Total Cases",
    y = "Observed Total Cases (2001–2023)"
  ) +
  theme_minimal()






#########################################
# Load libraries
#########################################
library(tidyverse)
library(janitor)
library(car)

#########################################
# --- Read and clean SVI data (3 vars) ---
#########################################
svi_all <- read_csv("data_raw/california_svi_2020.csv") %>%
  clean_names()

svi_summary <- svi_all %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    ep_pov150 = mean(ep_pov150, na.rm = TRUE),
    ep_disabl = mean(ep_disabl, na.rm = TRUE),
    ep_crowd = mean(ep_crowd, na.rm = TRUE),
    .groups = "drop"
  )

#########################################
# --- Read and clean PM2.5 data ---
#########################################
pm_raw <- read_csv("HDPulse_data_export.csv", col_names = "X1")
pm_trim <- pm_raw[-c(1:6), , drop = FALSE]

pm_clean <- pm_trim %>%
  separate(X1, into = c("county", "fips", "pm25"), sep = ",") %>%
  mutate(
    county = str_remove(county, " County"),
    county = str_to_upper(county),
    county = str_trim(county),
    pm25 = as.numeric(pm25)
  ) %>%
  filter(!is.na(pm25)) %>%
  dplyr::select(county, pm25)

#########################################
# --- Merge SVI and PM2.5 data ---
#########################################
svi_pm <- left_join(svi_summary, pm_clean, by = "county") %>%
  drop_na()

#########################################
# --- Read and clean Valley Fever data ---
#########################################
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),
    county = str_remove(county, " COUNTY"),
    county = str_trim(county),
    cases = as.numeric(cases_raw),
    county = case_when(
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  filter(total_cases > 44)

#########################################
# --- Merge all data together ---
#########################################
vf_final <- left_join(vf_clean, svi_pm, by = "county") %>%
  drop_na()

#########################################
# --- Build regression model (5 vars) ---
#########################################
final_formula <- "total_cases ~ ep_pov150 + ep_disabl + ep_crowd + pm25"
model_final <- lm(as.formula(final_formula), data = vf_final)

# Print summary
summary(model_final)

# Check VIFs
vif_values <- vif(model_final)
cat("VIF for each variable:\n")
print(vif_values)

#########################################
# --- Residual diagnostic plots ---
#########################################
par(mfrow = c(2, 2))
plot(model_final)
par(mfrow = c(1, 1))

#########################################
# --- Null model and AIC comparison ---
#########################################
model_null <- lm(total_cases ~ 1, data = vf_final)
aic_values <- AIC(model_final, model_null)
print(aic_values)

#########################################
# --- Predicted vs actual correlation ---
#########################################
vf_final <- vf_final %>%
  mutate(predicted_cases = predict(model_final, newdata = vf_final))

cor_test <- cor.test(vf_final$predicted_cases, vf_final$total_cases)
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)

# Print correlation explicitly
cat("Correlation (predicted vs. actual):", r_value, "\n")
cat("Correlation p-value:", p_value, "\n")

annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

#########################################
# --- Plot predicted vs actual ---
#########################################
ggplot(vf_final, aes(x = predicted_cases, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate("text",
           x = min(vf_final$predicted_cases, na.rm = TRUE) + 0.02,
           y = max(vf_final$total_cases, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4, color = "black") +
  labs(
    title = "Predicted vs. Actual Valley Fever Cases",
    x = "Predicted Cases (Fitted)",
    y = "Observed Cases (2001–2023)"
  ) +
  theme_minimal()















#########################################
# Load libraries
#########################################
library(tidyverse)
library(janitor)
library(mgcv)

#########################################
# --- Read and clean SVI data ---
#########################################
svi_all <- read_csv("data_raw/california_svi_2020.csv") %>%
  clean_names()

svi_summary <- svi_all %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    rpl_themes = mean(rpl_themes, na.rm = TRUE),
    e_totpop = mean(e_totpop, na.rm = TRUE),  # ✅ Add total population column
    .groups = "drop"
  )

#########################################
# --- Read and clean Valley Fever data ---
#########################################
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),
    county = str_remove(county, " COUNTY"),
    county = str_trim(county),
    cases = as.numeric(cases_raw),
    county = case_when(
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  filter(total_cases > 44)

#########################################
# --- Merge VF data and SVI summary ---
#########################################
vf_svi <- left_join(vf_clean, svi_summary, by = "county") %>%
  filter(!is.na(rpl_themes) & !is.na(e_totpop))

#########################################
# --- GAM model with offset ---
#########################################
model_gam <- gam(total_cases ~ s(rpl_themes) + offset(log(e_totpop)), family = nb(), data = vf_svi)

# Summary
summary(model_gam)

#########################################
# --- Residual plot ---
#########################################
res <- residuals(model_gam, type = "pearson")
plot(fitted(model_gam), res,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     main = "Residual plot (GAM with Offset)")
abline(h = 0, col = "red")

#########################################
# --- Predicted vs actual plot ---
#########################################
vf_svi <- vf_svi %>%
  mutate(predicted_cases_gam = predict(model_gam, type = "response"))

cor_test <- cor.test(vf_svi$predicted_cases_gam, vf_svi$total_cases)
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

ggplot(vf_svi, aes(x = predicted_cases_gam, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  annotate("text",
           x = min(vf_svi$predicted_cases_gam, na.rm = TRUE),
           y = max(vf_svi$total_cases, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4, color = "black") +
  labs(
    title = "Predicted vs. Actual Total Cases (GAM with Offset)",
    x = "Predicted Total Cases",
    y = "Observed Total Cases (2001–2023)"
  ) +
  theme_minimal()

#########################################
# --- Smooth function plot ---
#########################################
plot(model_gam, pages = 1, shade = TRUE)







#########################################
# Load libraries
#########################################
library(tidyverse)
library(janitor)
library(mgcv)

#########################################
# --- Read and clean SVI data ---
#########################################
svi_all <- read_csv("data_raw/california_svi_2020.csv") %>%
  clean_names()

svi_summary <- svi_all %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    pov = mean(ep_pov150, na.rm = TRUE),
    disabl = mean(ep_disabl, na.rm = TRUE),
    crowd = mean(ep_crowd, na.rm = TRUE),
    .groups = "drop"
  )

#########################################
# --- Read and clean Valley Fever data ---
#########################################
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),
    county = str_remove(county, " COUNTY"),
    county = str_trim(county),
    cases = as.numeric(cases_raw),
    county = case_when(
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  filter(total_cases > 44)

#########################################
# --- Merge data ---
#########################################
vf_svi <- left_join(vf_clean, svi_summary, by = "county") %>%
  filter(!is.na(pov) & !is.na(disabl) & !is.na(crowd))

#########################################
# --- GAM model ---
#########################################
model_gam_multi <- gam(total_cases ~ s(pov) + s(disabl) + s(crowd), family = nb(), data = vf_svi)

# Summary
summary(model_gam_multi)

#########################################
# --- Residual plot ---
#########################################
res_nb <- residuals(model_gam_multi, type = "pearson")
plot(fitted(model_gam_multi), res_nb,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     main = "Residual plot (GAM SVI 3 vars)")
abline(h = 0, col = "red")

#########################################
# --- Predicted vs actual plot ---
#########################################
vf_svi <- vf_svi %>%
  mutate(predicted_cases_gam = predict(model_gam_multi, type = "response"))

cor_test <- cor.test(vf_svi$predicted_cases_gam, vf_svi$total_cases)
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

ggplot(vf_svi, aes(x = predicted_cases_gam, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  annotate("text",
           x = min(vf_svi$predicted_cases_gam, na.rm = TRUE),
           y = max(vf_svi$total_cases, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4, color = "black") +
  labs(
    title = "Predicted vs. Actual Total Cases (GAM SVI 3 vars)",
    x = "Predicted Total Cases",
    y = "Observed Total Cases (2001–2023)"
  ) +
  theme_minimal()

#########################################
# --- Smooth function plots ---
#########################################
plot(model_gam_multi, pages = 1, shade = TRUE)






#########################################
# Load libraries
#########################################
library(tidyverse)
library(janitor)
library(mgcv)

#########################################
# --- Read and clean SVI data ---
#########################################
svi_all <- read_csv("data_raw/california_svi_2020.csv") %>%
  clean_names()

svi_summary <- svi_all %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    ep_pov150 = mean(ep_pov150, na.rm = TRUE),
    ep_disabl = mean(ep_disabl, na.rm = TRUE),
    ep_crowd = mean(ep_crowd, na.rm = TRUE),
    .groups = "drop"
  )

#########################################
# --- Read and clean PM2.5 data ---
#########################################
pm_raw <- read_csv("HDPulse_data_export.csv", col_names = "X1")
pm_trim <- pm_raw[-c(1:6), , drop = FALSE]

pm_clean <- pm_trim %>%
  separate(X1, into = c("county", "fips", "pm25"), sep = ",") %>%
  mutate(
    county = str_remove(county, " County"),
    county = str_to_upper(county),
    county = str_trim(county),
    pm25 = as.numeric(pm25)
  ) %>%
  filter(!is.na(pm25)) %>%
  dplyr::select(county, pm25)

#########################################
# --- Merge SVI and PM2.5 data ---
#########################################
svi_pm <- left_join(svi_summary, pm_clean, by = "county") %>%
  drop_na()

#########################################
# --- Read and clean Valley Fever data ---
#########################################
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),
    county = str_remove(county, " COUNTY"),
    county = str_trim(county),
    cases = as.numeric(cases_raw),
    county = case_when(
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  filter(total_cases > 44)

#########################################
# --- Merge final data ---
#########################################
vf_final <- left_join(vf_clean, svi_pm, by = "county") %>%
  drop_na()

#########################################
# --- GAM model ---
#########################################
model_gam_final <- gam(total_cases ~ s(ep_pov150) + s(ep_disabl) + s(ep_crowd) + s(pm25), family = nb(), data = vf_final)

# Summary
summary(model_gam_final)

#########################################
# --- Residual plot ---
#########################################
res_nb <- residuals(model_gam_final, type = "pearson")
plot(fitted(model_gam_final), res_nb,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     main = "Residual plot (GAM SVI + PM2.5)")
abline(h = 0, col = "red")

#########################################
# --- Predicted vs actual plot ---
#########################################
vf_final <- vf_final %>%
  mutate(predicted_cases_gam = predict(model_gam_final, type = "response"))

cor_test <- cor.test(vf_final$predicted_cases_gam, vf_final$total_cases)
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

ggplot(vf_final, aes(x = predicted_cases_gam, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  annotate("text",
           x = min(vf_final$predicted_cases_gam, na.rm = TRUE),
           y = max(vf_final$total_cases, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4, color = "black") +
  labs(
    title = "Predicted vs. Actual Total Cases (GAM SVI + PM2.5)",
    x = "Predicted Total Cases",
    y = "Observed Total Cases (2001–2023)"
  ) +
  theme_minimal()

#########################################
# --- Smooth function plots ---
#########################################
plot(model_gam_final, pages = 1, shade = TRUE)

