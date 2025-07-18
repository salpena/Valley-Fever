
#########################################
# --- Install required packages (only once) ---
#########################################
# Uncomment and run this block once if packages are not installed
# install.packages(c("tidyverse", "janitor", "car", "mgcv", "leaps", "corrplot"))

#########################################
# --- Load libraries ---
#########################################
library(tidyverse)   # Data manipulation and plotting
library(janitor)     # Clean column names
library(car)         # Variance Inflation Factor (VIF)
library(mgcv)        # Generalized Additive Models (GAM)
library(leaps)       # Best subset regression
library(corrplot)    # Correlation matrix visualization

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

# Count number of counties with >44 total cases
num_counties <- vf_clean %>% nrow()

# Print the number
cat("Number of counties with >44 total cases:", num_counties, "\n")


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
# Clean and prepare SVI data.
# Calculate correlations among all 16 SVI variables and plot a heatmap to check for strong linear relationships.
# Calculate Variance Inflation Factors (VIF) to detect multicollinearity among predictors — high VIF means a variable is highly redundant.
# Visualize VIFs using a bar plot with threshold lines at 5 and 10 for easy interpretation.



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

# This code first cleans and summarizes SVI data by county.
# It then cleans Valley Fever case data and merges with SVI data.
# We use best subset regression to try every possible combination of SVI predictors (up to 16).
# For each combination, it checks if all VIF values are below 5 to avoid multicollinearity.
# Successful models (with VIF < 5) and their adjusted R² values are printed to help select the most stable models.


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
  filter(total_cases > 44)  #Keep only counties with more than 44 cases

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
# --- Read and clean SVI data (16 variables) ---
#########################################
# Load raw SVI dataset and standardize column names to snake_case
svi_all <- read_csv("data_raw/california_svi_2020.csv") %>%
  clean_names()

# Calculate average of each SVI component across tracts by county
svi_summary <- svi_all %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    ep_pov150 = mean(ep_pov150, na.rm = TRUE),   # % below 150% of poverty line
    ep_unemp = mean(ep_unemp, na.rm = TRUE),     # % unemployed
    ep_hburd = mean(ep_hburd, na.rm = TRUE),     # % with housing cost burden
    ep_nohsdp = mean(ep_nohsdp, na.rm = TRUE),   # % without high school diploma
    ep_uninsur = mean(ep_uninsur, na.rm = TRUE), # % uninsured
    ep_age65 = mean(ep_age65, na.rm = TRUE),     # % over age 65
    ep_age17 = mean(ep_age17, na.rm = TRUE),     # % under age 17
    ep_disabl = mean(ep_disabl, na.rm = TRUE),   # % disabled
    ep_sngpnt = mean(ep_sngpnt, na.rm = TRUE),   # % single-parent households
    ep_limeng = mean(ep_limeng, na.rm = TRUE),   # % with limited English
    ep_minrty = mean(ep_minrty, na.rm = TRUE),   # % minority
    ep_munit = mean(ep_munit, na.rm = TRUE),     # % in multi-unit housing
    ep_mobile = mean(ep_mobile, na.rm = TRUE),   # % in mobile homes
    ep_crowd = mean(ep_crowd, na.rm = TRUE),     # % living in crowded housing
    ep_noveh = mean(ep_noveh, na.rm = TRUE),      # % without vehicles
    ep_groupq = mean(ep_groupq, na.rm = TRUE),    # % in group quarters
    .groups = "drop"
  )

#########################################
# --- Read and clean PM2.5 data ---
#########################################
# Read raw PM2.5 CSV file (downloaded from HDPulse)
pm_raw <- read_csv("HDPulse_data_export.csv", col_names = "X1")

# Remove non-data header and national average rows
pm_trim <- pm_raw[-c(1:6), , drop = FALSE]

# Split single-column string into multiple variables
pm_clean <- pm_trim %>%
  separate(X1, into = c("county", "fips", "pm25"), sep = ",") %>%
  mutate(
    county = str_remove(county, " County"),      # Remove "County" text
    county = str_to_upper(county),               # Match format with other datasets
    county = str_trim(county),
    pm25 = as.numeric(pm25)                      # Convert to numeric
  ) %>%
  filter(!is.na(pm25)) %>%
  dplyr::select(county, pm25)                    # Keep only needed columns

#########################################
# --- Merge SVI and PM2.5 data ---
#########################################
# Combine SVI summary data with air pollution data by county
svi_pm <- left_join(svi_summary, pm_clean, by = "county") %>%
  drop_na()

#########################################
# --- Read and clean Valley Fever data ---
#########################################
# Load Valley Fever case data and assign column names
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

# Clean VF data and aggregate total cases per county
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
# Join VF case counts with SVI and PM2.5 variables
vf_full <- left_join(vf_clean, svi_pm, by = "county") %>%
  drop_na()

#########################################
# --- Run best subset regression (all 17 variables) ---
#########################################
# Create formula: total_cases ~ all 17 predictors
predictor_vars <- names(vf_full)[!(names(vf_full) %in% c("county", "total_cases"))]
formula_text <- paste("total_cases ~", paste(predictor_vars, collapse = " + "))

# Run best subset regression using regsubsets()
best_subset <- regsubsets(as.formula(formula_text), data = vf_full, nvmax = length(predictor_vars), really.big = TRUE)
summary_best <- summary(best_subset)

#########################################
# --- Check combinations with VIF < 5 ---
#########################################
# Identify variable combos with acceptable multicollinearity (VIF < 5)
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
# Display model combos that passed the VIF < 5 filter
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
# --- Read and clean Valley Fever data ---
#########################################
# Read raw Valley Fever case data, skipping headers
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

# Clean and preprocess the VF data
vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),                  # Standardize county names to uppercase
    county = str_remove(county, " COUNTY"),         # Remove " COUNTY" suffix
    county = str_trim(county),                      # Remove leading/trailing spaces
    cases = as.numeric(cases_raw),                  # Convert case counts to numeric
    county = case_when(                             # Reassign some cities to main counties
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%  # Remove aggregate/unknown rows
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  filter(total_cases > 44)  # Keep counties with higher case burden

#########################################
# --- Merge all data together ---
#########################################
# Join cleaned VF data with selected SVI + PM2.5 predictors
vf_final <- left_join(vf_clean, svi_pm, by = "county") %>%
  drop_na()

#########################################
# --- Build regression model using 4 vars ---
#########################################
# Build multiple linear regression using 4 selected predictors
final_formula <- "total_cases ~ ep_pov150 + ep_disabl + ep_crowd + pm25"
model_final <- lm(as.formula(final_formula), data = vf_final)

# View regression coefficients, significance levels, and model fit
summary(model_final)

#########################################
# --- Check VIFs for multicollinearity ---
#########################################
# Check for multicollinearity among predictors using Variance Inflation Factor
vif_values <- vif(model_final)
cat("VIF for each variable:\n")
print(vif_values)

#########################################
# --- Add predicted values ---
#########################################
# Add predicted case values to the dataset
vf_final <- vf_final %>%
  mutate(predicted_cases = predict(model_final, newdata = vf_final))

#########################################
# --- Correlation test ---
#########################################
# Correlation between predicted and actual total cases
cor_test <- cor.test(vf_final$predicted_cases, vf_final$total_cases)
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)

cat("Correlation between predicted and actual cases (r):", r_value, "\n")
cat("p-value:", p_value, "\n")

# Format annotation for plot
annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

#########################################
# --- Plot predicted vs actual cases ---
#########################################
ggplot(vf_final, aes(x = predicted_cases, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Add regression line
  annotate("text",
           x = min(vf_final$predicted_cases, na.rm = TRUE) + 0.02,
           y = max(vf_final$total_cases, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4, color = "black") +
  labs(
    title = "Predicted vs. Actual Valley Fever Cases (4-Variable Model)",
    x = "Predicted Cases",
    y = "Observed Cases (2001–2023)"
  ) +
  theme_minimal()

#########################################
# --- Residual diagnostic plots ---
#########################################
# Generate default residual diagnostic plots (QQ, residuals vs fitted, etc.)
par(mfrow = c(2, 2))  # 2x2 panel
plot(model_final)
par(mfrow = c(1, 1))  # Reset layout

#########################################
# --- Null model and AIC comparison ---
#########################################
# Compare model fit with null model using AIC
model_null <- lm(total_cases ~ 1, data = vf_final)
AIC(model_final, model_null)

#########################################
# --- Read and clean SVI data ---
#########################################
# Clean county-level SVI overall theme score and population data
svi_all <- read_csv("data_raw/california_svi_2020.csv") %>%
  clean_names()

svi_summary <- svi_all %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    rpl_themes = mean(rpl_themes, na.rm = TRUE),     # SVI composite score
    e_totpop = mean(e_totpop, na.rm = TRUE),         # Total population
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
# Merge VF and SVI data for GAM model
vf_svi <- left_join(vf_clean, svi_summary, by = "county") %>%
  filter(!is.na(rpl_themes) & !is.na(e_totpop))

#########################################
# --- GAM model with offset ---
#########################################
# Generalized Additive Model using negative binomial distribution
# Include offset for total population to control for exposure
model_gam <- gam(total_cases ~ s(rpl_themes) + offset(log(e_totpop)), family = nb(), data = vf_svi)

# View GAM summary: smooth term fit, significance, deviance explained
summary(model_gam)

#########################################
# --- Residual plot ---
#########################################
# Plot Pearson residuals to check model fit visually
res <- residuals(model_gam, type = "pearson")
plot(fitted(model_gam), res,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     main = "Residual plot (GAM with Offset)")
abline(h = 0, col = "red")

#########################################
# --- Predicted vs actual plot ---
#########################################
# Generate predicted values from GAM and compare with observed
vf_svi <- vf_svi %>%
  mutate(predicted_cases_gam = predict(model_gam, type = "response"))

cor_test <- cor.test(vf_svi$predicted_cases_gam, vf_svi$total_cases)
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

# Plot predicted vs observed cases (GAM)
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
# Visualize smooth term (spline for SVI score)
plot(model_gam, pages = 1, shade = TRUE)

#########################################
# --- Compute and print AIC ---
#########################################
# Model selection using AIC: lower is better
model_aic <- AIC(model_gam)
cat("AIC for GAM with offset:", model_aic, "\n")



#########################################
# --- Read and clean SVI data ---
#########################################

# Load SVI dataset and clean column names
svi_all <- read_csv("data_raw/california_svi_2020.csv") %>%
  clean_names()

# Summarize SVI variables by county: calculate mean of key indicators
svi_summary <- svi_all %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    pov = mean(ep_pov150, na.rm = TRUE),       # Poverty rate
    disabl = mean(ep_disabl, na.rm = TRUE),    # Disability rate
    crowd = mean(ep_crowd, na.rm = TRUE),      # Crowded housing
    .groups = "drop"
  )

#########################################
# --- Read and clean Valley Fever data ---
#########################################

# Load Valley Fever cases dataset, skip metadata rows, assign column names
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

# Clean county names, fix special jurisdictions, and summarize case totals
vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),
    county = str_remove(county, " COUNTY"),
    county = str_trim(county),
    cases = as.numeric(cases_raw),
    county = case_when(     # Combine independent cities into respective counties
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%  # Remove total rows or asterisks
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  filter(total_cases > 44)  # Filter to counties with more than 44 cases

#########################################
# --- Merge data ---
#########################################

# Join VF and SVI data on county name and remove any rows with missing SVI info
vf_svi <- left_join(vf_clean, svi_summary, by = "county") %>%
  filter(!is.na(pov) & !is.na(disabl) & !is.na(crowd))

#########################################
# --- GAM model (SVI only) ---
#########################################

# Fit a Generalized Additive Model (GAM) with a Negative Binomial family
model_gam_multi <- gam(total_cases ~ s(pov) + s(disabl) + s(crowd),
                       family = nb(), data = vf_svi)

# Print model summary (estimates, significance, smoothness)
summary(model_gam_multi)

#########################################
# --- Residual plot ---
#########################################

# Plot Pearson residuals to check model fit
res_nb <- residuals(model_gam_multi, type = "pearson")
plot(fitted(model_gam_multi), res_nb,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     main = "Residual plot (GAM SVI 3 vars)")
abline(h = 0, col = "red")  # Reference line at 0

#########################################
# --- Predicted vs actual plot ---
#########################################

# Generate predicted case counts and compute correlation with observed cases
vf_svi <- vf_svi %>%
  mutate(predicted_cases_gam = predict(model_gam_multi, type = "response"))

cor_test <- cor.test(vf_svi$predicted_cases_gam, vf_svi$total_cases)
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

# Plot predicted vs. actual case totals with correlation annotation
ggplot(vf_svi, aes(x = predicted_cases_gam, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = min(vf_svi$predicted_cases_gam, na.rm = TRUE),
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

# Plot smooth terms (nonlinear relationships for each predictor)
plot(model_gam_multi, pages = 1, shade = TRUE)

#########################################
# --- Compute and print AIC ---
#########################################

# AIC for model selection (lower = better fit)
model_aic <- AIC(model_gam_multi)
cat("AIC for GAM with 3 SVI variables:", model_aic, "\n")

##########################################################
# --- Repeat with additional environmental data: PM2.5 ---
##########################################################

# Load and summarize SVI data again with original column names
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

# Read and clean PM2.5 data from a CSV where all data is in a single column
pm_raw <- read_csv("HDPulse_data_export.csv", col_names = "X1")
pm_trim <- pm_raw[-c(1:6), , drop = FALSE]  # Remove header rows

# Split column into multiple fields and clean up
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

# Merge PM2.5 data with SVI summary
svi_pm <- left_join(svi_summary, pm_clean, by = "county") %>%
  drop_na()

# Reload and clean Valley Fever data again for merge
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

# Merge all data (VF cases + SVI + PM2.5)
vf_final <- left_join(vf_clean, svi_pm, by = "county") %>%
  drop_na()

#########################################
# --- GAM model (SVI + PM2.5) ---
#########################################

# Fit new GAM including air pollution (PM2.5)
model_gam_final <- gam(total_cases ~ s(ep_pov150) + s(ep_disabl) + s(ep_crowd) + s(pm25),
                       family = nb(), data = vf_final)

# View summary
summary(model_gam_final)

# Residual plot
res_nb <- residuals(model_gam_final, type = "pearson")
plot(fitted(model_gam_final), res_nb,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     main = "Residual plot (GAM SVI + PM2.5)")
abline(h = 0, col = "red")

# Predicted vs. actual
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

# Smooth plots for nonlinear terms
plot(model_gam_final, pages = 1, shade = TRUE)

# AIC for model comparison
model_aic <- AIC(model_gam_final)
cat("AIC for GAM with SVI + PM2.5 variables:", model_aic, "\n")



#########################################
# --- Read individual yearly temperature files ---
#########################################

# Define folder containing the yearly temperature CSVs
folder_path <- "Temp"

# Read each year's temperature CSV individually into separate objects (manual method)
# These objects are not used in the cleaned analysis but kept for reference/testing
# Read each CSV individually and assign to separate objects
data_2001 <- read_csv(file.path(folder_path, "2001.csv"))
data_2002 <- read_csv(file.path(folder_path, "2002.csv"))
data_2003 <- read_csv(file.path(folder_path, "2003.csv"))
data_2004 <- read_csv(file.path(folder_path, "2004.csv"))
data_2005 <- read_csv(file.path(folder_path, "2005.csv"))
data_2006 <- read_csv(file.path(folder_path, "2006.csv"))
data_2007 <- read_csv(file.path(folder_path, "2007.csv"))
data_2008 <- read_csv(file.path(folder_path, "2008.csv"))
data_2009 <- read_csv(file.path(folder_path, "2009.csv"))
data_2010 <- read_csv(file.path(folder_path, "2010.csv"))
data_2011 <- read_csv(file.path(folder_path, "2011.csv"))
data_2012 <- read_csv(file.path(folder_path, "2012.csv"))
data_2013 <- read_csv(file.path(folder_path, "2013.csv"))
data_2014 <- read_csv(file.path(folder_path, "2014.csv"))
data_2015 <- read_csv(file.path(folder_path, "2015.csv"))
data_2016 <- read_csv(file.path(folder_path, "2016.csv"))
data_2017 <- read_csv(file.path(folder_path, "2017.csv"))
data_2018 <- read_csv(file.path(folder_path, "2018.csv"))
data_2019 <- read_csv(file.path(folder_path, "2019.csv"))
data_2020 <- read_csv(file.path(folder_path, "2020.csv"))
data_2021 <- read_csv(file.path(folder_path, "2021.csv"))
data_2022 <- read_csv(file.path(folder_path, "2022.csv"))
data_2023 <- read_csv(file.path(folder_path, "2023.csv"))

#########################################
# --- Loop through temperature files and clean ---
#########################################

library(tidyverse)

# Define the vector of years to automate file reading
years <- 2001:2023
all_data <- list()  # Initialize an empty list to store cleaned data

# Loop over each year and clean temperature data
for (yr in years) {
  file_path <- paste0("Temp/", yr, ".csv")
  temp_raw <- read_csv(file_path, skip = 4)  # Skip metadata/header rows
  
  temp_clean <- temp_raw %>%
    dplyr::select(Name, Value) %>%  # Keep only relevant columns
    mutate(year = yr)               # Add year column
  
  all_data[[as.character(yr)]] <- temp_clean  # Store in list
}

# Combine all yearly datasets into a single dataframe
combined_data <- bind_rows(all_data)

# Convert temperature values to numeric (they may be read as character)
combined_data <- combined_data %>%
  mutate(Value = as.numeric(Value))

# Compute average temperature per county over 2001–2023
county_avg <- combined_data %>%
  group_by(Name) %>%
  summarise(
    total_value = sum(Value, na.rm = TRUE),
    avg_value = total_value / 22  # Average over 22 years
  ) %>%
  arrange(desc(avg_value))  # Optional: sort by average temp

# View results
print(county_avg)

#########################################
# --- Read and clean Valley Fever data ---
#########################################

# Read VF case data, clean column names, and fix city/county entries
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),                    # Standardize capitalization
    county = str_remove(county, " COUNTY"),           # Remove "County" label
    county = str_trim(county),
    cases = as.numeric(cases_raw),
    county = case_when(                               # Combine sub-county jurisdictions
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  filter(total_cases > 44)  # Keep counties with meaningful case totals

#########################################
# --- Clean and merge temperature data ---
#########################################

# Standardize county names in temperature data
county_avg <- county_avg %>%
  mutate(Name = str_to_upper(Name),
         Name = str_remove(Name, " COUNTY"),
         Name = str_trim(Name))

# Merge VF and average temperature data for valid counties
merged_vf_temp <- vf_clean %>%
  left_join(county_avg, by = c("county" = "Name")) %>%
  drop_na(avg_value)

# View final dataset
print(merged_vf_temp)

#########################################
# --- Pearson correlation analysis ---
#########################################

# Check strength and significance of correlation
cor_test <- cor.test(merged_vf_temp$total_cases, merged_vf_temp$avg_value)

# Display correlation results
cat("Pearson correlation coefficient (r):", round(cor_test$estimate, 3), "\n")
cat("p-value:", signif(cor_test$p.value, 3), "\n")

# Prepare annotation for scatter plot
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

# Plot VF cases vs. average temperature
ggplot(merged_vf_temp, aes(x = avg_value, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate("text",
           x = min(merged_vf_temp$avg_value, na.rm = TRUE) + 0.2,
           y = max(merged_vf_temp$total_cases, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4, color = "black") +
  labs(
    title = "Total Valley Fever Cases vs. Average Temperature (2001–2023)",
    x = "Average Temperature (°F)",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()

#########################################
# --- Prepare full merged dataset ---
#########################################

# Re-clean and merge SVI and PM2.5 data (same as prior script)
# See earlier documentation for these steps...

# Add temperature data to SVI + PM2.5 + VF merged dataset
vf_final <- left_join(vf_final,
                      merged_vf_temp %>% dplyr::select(county, avg_value),
                      by = "county") %>%
  rename(avg_temp = avg_value) %>%
  drop_na()

#########################################
# --- Build 5-variable linear regression model ---
#########################################

# Define formula with 3 SVI vars + PM2.5 + temperature
final_formula <- "total_cases ~ ep_pov150 + ep_disabl + ep_crowd + pm25 + avg_temp"
model_final <- lm(as.formula(final_formula), data = vf_final)

# Print model summary
summary(model_final)

#########################################
# --- Check for multicollinearity ---
#########################################

# Calculate Variance Inflation Factors (VIF) to check collinearity
vif_values <- car::vif(model_final)
cat("VIF for each variable:\n")
print(vif_values)

#########################################
# --- Predict total VF cases ---
#########################################

# Add predicted values from linear model to dataset
vf_final <- vf_final %>%
  mutate(predicted_cases = predict(model_final, newdata = vf_final))

#########################################
# --- Correlation test (predicted vs. actual) ---
#########################################

# Assess model predictive strength
cor_test <- cor.test(vf_final$predicted_cases, vf_final$total_cases)
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
cat("Correlation between predicted and actual cases (r):", r_value, "\n")
cat("p-value:", p_value, "\n")

annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

#########################################
# --- Plot predicted vs. actual cases ---
#########################################

# Scatter plot with trend line and correlation annotation
ggplot(vf_final, aes(x = predicted_cases, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate("text",
           x = min(vf_final$predicted_cases, na.rm = TRUE) + 0.02,
           y = max(vf_final$total_cases, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4, color = "black") +
  labs(
    title = "Predicted vs. Actual Valley Fever Cases (5-Variable Model)",
    x = "Predicted Cases",
    y = "Observed Cases (2001–2023)"
  ) +
  theme_minimal()

#########################################
# --- Residual diagnostics ---
#########################################

# Plot diagnostic checks for linear model assumptions
par(mfrow = c(2, 2))
plot(model_final)
par(mfrow = c(1, 1))  # Reset to default layout

#########################################
# --- AIC model selection ---
#########################################

# Print AIC value for model comparison (lower is better)
model_aic <- AIC(model_final)
cat("AIC for 5-variable linear model:", model_aic, "\n")

#########################################
# --- Fit GAM with all 5 variables ---
#########################################

# Build a flexible model using nonlinear smooths (GAM)
model_gam_final <- gam(total_cases ~ 
                         s(ep_pov150) + 
                         s(ep_disabl) + 
                         s(ep_crowd) + 
                         s(pm25) + 
                         s(avg_temp),
                       family = nb(), data = vf_final)

# Model summary
summary(model_gam_final)

# Residual plot for GAM
res_nb <- residuals(model_gam_final, type = "pearson")
plot(fitted(model_gam_final), res_nb,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     main = "Residual plot (GAM SVI + PM2.5 + Temp)")
abline(h = 0, col = "red")

# Predicted vs. actual for GAM
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
    title = "Predicted vs. Actual Total Cases (GAM SVI + PM2.5 + Temp)",
    x = "Predicted Total Cases",
    y = "Observed Total Cases (2001–2023)"
  ) +
  theme_minimal()

# Plot smooth terms to visualize nonlinear effects
plot(model_gam_final, pages = 1, shade = TRUE)

# AIC for GAM model
model_aic <- AIC(model_gam_final)
cat("AIC for GAM with SVI + PM2.5 + Temp variables:", model_aic, "\n")




# Folder path where your files are stored
folder_path <- "Water"

# Read each CSV individually and assign to separate objects
data_p2001 <- read_csv(file.path(folder_path, "p2001.csv"), skip = 4)
data_p2002 <- read_csv(file.path(folder_path, "p2002.csv"), skip = 4)
data_p2003 <- read_csv(file.path(folder_path, "p2003.csv"), skip = 4)
data_p2004 <- read_csv(file.path(folder_path, "p2004.csv"), skip = 4)
data_p2005 <- read_csv(file.path(folder_path, "p2005.csv"), skip = 4)
data_p2006 <- read_csv(file.path(folder_path, "p2006.csv"), skip = 4)
data_p2007 <- read_csv(file.path(folder_path, "p2007.csv"), skip = 4)
data_p2008 <- read_csv(file.path(folder_path, "p2008.csv"), skip = 4)
data_p2009 <- read_csv(file.path(folder_path, "p2009.csv"), skip = 4)
data_p2010 <- read_csv(file.path(folder_path, "p2010.csv"), skip = 4)
data_p2011 <- read_csv(file.path(folder_path, "p2011.csv"), skip = 4)
data_p2012 <- read_csv(file.path(folder_path, "p2012.csv"), skip = 4)
data_p2013 <- read_csv(file.path(folder_path, "p2013.csv"), skip = 4)
data_p2014 <- read_csv(file.path(folder_path, "p2014.csv"), skip = 4)
data_p2015 <- read_csv(file.path(folder_path, "p2015.csv"), skip = 4)
data_p2016 <- read_csv(file.path(folder_path, "p2016.csv"), skip = 4)
data_p2017 <- read_csv(file.path(folder_path, "p2017.csv"), skip = 4)
data_p2018 <- read_csv(file.path(folder_path, "p2018.csv"), skip = 4)
data_p2019 <- read_csv(file.path(folder_path, "p2019.csv"), skip = 4)
data_p2020 <- read_csv(file.path(folder_path, "p2020.csv"), skip = 4)
data_p2021 <- read_csv(file.path(folder_path, "p2021.csv"), skip = 4)
data_p2022 <- read_csv(file.path(folder_path, "p2022.csv"), skip = 4)
data_p2023 <- read_csv(file.path(folder_path, "p2023.csv"), skip = 4)


#########################################
# --- Loop through precipitation files and clean ---
#########################################

# Define vector of years to loop over
years_p <- 2001:2023
all_data_p <- list()  # Initialize empty list to store cleaned annual data

# Loop over each year’s precipitation CSV
for (yr in years_p) {
  file_path_p <- paste0("Water/p", yr, ".csv")  # Build file path
  precip_raw_p <- read_csv(file_path_p, skip = 4)  # Skip metadata rows
  
  precip_clean_p <- precip_raw_p %>%
    dplyr::select(Name, Value) %>%  # Keep only county and precipitation value
    mutate(year = yr)               # Add year column
  
  all_data_p[[as.character(yr)]] <- precip_clean_p  # Save in list
}

# Combine all yearly datasets into one dataframe
combined_data_p <- bind_rows(all_data_p)

# Convert Value column to numeric (in case it's read as character)
combined_data_p <- combined_data_p %>%
  mutate(Value = as.numeric(Value))

#########################################
# --- Calculate average precipitation per county ---
#########################################

# Group by county and calculate total and average precipitation
county_avg_p <- combined_data_p %>%
  group_by(Name) %>%
  summarise(
    total_value_p = sum(Value, na.rm = TRUE),        # Total over 22 years
    avg_value_p = total_value_p / 22                 # Annual average
  ) %>%
  arrange(desc(avg_value_p))  # Optional: sort by average value

# View summary
print(county_avg_p)

#########################################
# --- Read and clean Valley Fever data ---
#########################################

vf_cases_p <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                       col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean_p <- vf_cases_p %>%
  mutate(
    county = str_to_upper(county),                      # Standardize capitalization
    county = str_remove(county, " COUNTY"),             # Remove "County" suffix
    county = str_trim(county),
    cases = as.numeric(cases_raw),
    county = case_when(                                 # Merge independent cities into counties
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%  # Remove summary rows
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  filter(total_cases > 44)  # Keep counties with meaningful data

#########################################
# --- Clean and merge precipitation data ---
#########################################

# Clean county names in precipitation dataset
county_avg_p <- county_avg_p %>%
  mutate(Name = str_to_upper(Name),
         Name = str_remove(Name, " COUNTY"),
         Name = str_trim(Name))

# Merge VF data with precipitation averages by county
merged_vf_precip_p <- vf_clean_p %>%
  left_join(county_avg_p, by = c("county" = "Name")) %>%
  drop_na(avg_value_p)  # Drop counties missing precipitation data

# View final merged data
print(merged_vf_precip_p)

#########################################
# --- Pearson correlation analysis ---
#########################################

# Test correlation between average precipitation and VF cases
cor_test_p <- cor.test(merged_vf_precip_p$total_cases, merged_vf_precip_p$avg_value_p)

# Print r and p values
cat("Pearson correlation coefficient (r):", round(cor_test_p$estimate, 3), "\n")
cat("p-value:", signif(cor_test_p$p.value, 3), "\n")

# Format annotation text for plotting
r_value_p <- round(cor_test_p$estimate, 3)
p_value_p <- signif(cor_test_p$p.value, 3)
annotation_text_p <- sprintf("r = %.3f\np = %.3f", r_value_p, p_value_p)

#########################################
# --- Plot VF cases vs. average precipitation ---
#########################################

ggplot(merged_vf_precip_p, aes(x = avg_value_p, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +                          # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "blue") +     # Linear trendline
  annotate("text",
           x = min(merged_vf_precip_p$avg_value_p, na.rm = TRUE) + 0.2,
           y = max(merged_vf_precip_p$total_cases, na.rm = TRUE) * 0.9,
           label = annotation_text_p, hjust = 0, size = 4, color = "black") +
  labs(
    title = "Total Valley Fever Cases vs. Average Precipitation (2001–2023)",
    x = "Average Precipitation (inches or mm)",              # Units depend on data source
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()


#########################################
# --- Read and clean Valley Fever data ---
#########################################

# Load Valley Fever dataset and clean county names
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(
    county = str_to_upper(county),                      # Standardize names
    county = str_remove(county, " COUNTY"),
    county = str_trim(county),
    cases = as.numeric(cases_raw),
    county = case_when(                                 # Merge sub-county cities
      county == "BERKELEY" ~ "ALAMEDA",
      county == "LONG BEACH" ~ "LOS ANGELES",
      county == "PASADENA" ~ "LOS ANGELES",
      TRUE ~ county
    )
  ) %>%
  filter(!is.na(cases) & !str_detect(county, "TOTAL|\\*")) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  filter(total_cases > 44)  # Keep counties with sufficient data

#########################################
# --- Merge SVI + PM2.5 + Temperature + Precipitation ---
#########################################

# Merge SVI + PM2.5 data
vf_final <- left_join(vf_clean, svi_pm, by = "county") %>%
  drop_na()

# Merge temperature
vf_final <- left_join(
  vf_final,
  merged_vf_temp %>% dplyr::select(county, avg_value),
  by = "county"
) %>%
  rename(avg_temp = avg_value) %>%
  drop_na()

# Merge precipitation
vf_final <- left_join(
  vf_final,
  merged_vf_precip_p %>% dplyr::select(county, avg_value_p),
  by = "county"
) %>%
  rename(avg_precip = avg_value_p) %>%
  drop_na()

#########################################
# --- Linear regression: 6-variable model ---
#########################################

# Fit a linear regression model with SVI + PM2.5 + Temp + Precip
final_formula <- "total_cases ~ ep_pov150 + ep_disabl + ep_crowd + pm25 + avg_temp + avg_precip"
model_final <- lm(as.formula(final_formula), data = vf_final)

# Print model coefficients and statistics
summary(model_final)

#########################################
# --- Multicollinearity check (VIF) ---
#########################################

# Calculate VIF to detect multicollinearity
vif_values <- car::vif(model_final)
cat("VIF for each variable:\n")
print(vif_values)

#########################################
# --- Predicted values and model correlation ---
#########################################

# Add predicted values to dataframe
vf_final <- vf_final %>%
  mutate(predicted_cases = predict(model_final, newdata = vf_final))

# Compute correlation between predicted and actual cases
cor_test <- cor.test(vf_final$predicted_cases, vf_final$total_cases)
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
cat("Correlation between predicted and actual cases (r):", r_value, "\n")
cat("p-value:", p_value, "\n")

# Format annotation text
annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

#########################################
# --- Plot: Predicted vs. Actual (Linear Model) ---
#########################################

ggplot(vf_final, aes(x = predicted_cases, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate("text",
           x = min(vf_final$predicted_cases, na.rm = TRUE) + 0.02,
           y = max(vf_final$total_cases, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4) +
  labs(
    title = "Predicted vs. Actual Valley Fever Cases (6-Variable Model)",
    x = "Predicted Cases",
    y = "Observed Cases (2001–2023)"
  ) +
  theme_minimal()

#########################################
# --- Residual diagnostic plots (Linear Model) ---
#########################################

par(mfrow = c(2, 2))
plot(model_final)
par(mfrow = c(1, 1))  # Reset plot layout

#########################################
# --- AIC comparison ---
#########################################

# AIC for model selection (lower = better fit)
model_aic <- AIC(model_final)
cat("AIC for 6-variable linear model:", model_aic, "\n")

#########################################
# --- GAM model: All 6 variables ---
#########################################

# Fit GAM with smoothing for each variable (allows nonlinear effects)
model_gam_6 <- gam(total_cases ~ 
                     s(ep_pov150) + s(ep_disabl) + s(ep_crowd) + 
                     s(pm25) + s(avg_temp) + s(avg_precip),
                   family = nb(), data = vf_final)

summary(model_gam_6)

# Residual plot
res_nb <- residuals(model_gam_6, type = "pearson")
plot(fitted(model_gam_6), res_nb,
     xlab = "Fitted values", ylab = "Pearson residuals",
     main = "Residual plot (GAM SVI + PM2.5 + Temp + Precip)")
abline(h = 0, col = "red")

# Predicted vs. Actual
vf_final <- vf_final %>%
  mutate(predicted_cases_gam_6 = predict(model_gam_6, type = "response"))

cor_test <- cor.test(vf_final$predicted_cases_gam_6, vf_final$total_cases)
r_value <- round(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
annotation_text <- sprintf("r = %.3f\np = %.3f", r_value, p_value)

ggplot(vf_final, aes(x = predicted_cases_gam_6, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  annotate("text",
           x = min(vf_final$predicted_cases_gam_6, na.rm = TRUE),
           y = max(vf_final$total_cases, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4) +
  labs(
    title = "Predicted vs. Actual Total Cases (GAM SVI + PM2.5 + Temp + Precip)",
    x = "Predicted Total Cases",
    y = "Observed Total Cases (2001–2023)"
  ) +
  theme_minimal()

# Plot smooth terms
plot(model_gam_6, pages = 1, shade = TRUE)

# AIC
model_aic <- AIC(model_gam_6)
cat("AIC for GAM with SVI + PM2.5 + Temp + Precip:", model_aic, "\n")

#########################################
# --- GAM with interaction term (Tensor Product) ---
#########################################

# Fit GAM with interaction between temp and precip using tensor product
model_gam_ti <- gam(
  total_cases ~ 
    s(ep_pov150) + s(ep_disabl) + s(ep_crowd) + s(pm25) +
    s(avg_temp) + s(avg_precip) +
    ti(avg_temp, avg_precip),  # Add interaction term
  family = nb(),
  data = vf_final
)

summary(model_gam_ti)
AIC(model_gam_ti)

# Residuals
res_ti <- residuals(model_gam_ti, type = "pearson")
plot(fitted(model_gam_ti), res_ti,
     xlab = "Fitted values", ylab = "Pearson residuals",
     main = "Residual Plot: GAM with ti(temp, precip)")
abline(h = 0, col = "red")

# Predicted values and correlation
vf_final <- vf_final %>%
  mutate(predicted_cases_gam_ti = predict(model_gam_ti, type = "response"))

cor_test_ti <- cor.test(vf_final$predicted_cases_gam_ti, vf_final$total_cases)
r_value_ti <- round(cor_test_ti$estimate, 3)
p_value_ti <- signif(cor_test_ti$p.value, 3)

# Plot
ggplot(vf_final, aes(x = predicted_cases_gam_ti, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  annotate("text",
           x = min(vf_final$predicted_cases_gam_ti, na.rm = TRUE),
           y = max(vf_final$total_cases, na.rm = TRUE) * 0.9,
           label = sprintf("r = %.3f\np = %.3f", r_value_ti, p_value_ti),
           hjust = 0, size = 4) +
  labs(
    title = "Predicted vs. Actual Total Cases (GAM + ti(temp, precip))",
    x = "Predicted Total Cases",
    y = "Observed Total Cases"
  ) +
  theme_minimal()

# Plot smooth terms and tensor interaction
plot(model_gam_ti, pages = 1, shade = TRUE)

#########################################
# --- Environmental-only GAM model ---
#########################################

# Fit GAM with just environmental predictors and interaction
model_gam_env_only <- gam(total_cases ~ 
                            s(pm25) + 
                            s(avg_temp) + 
                            s(avg_precip) + 
                            ti(avg_temp, avg_precip),
                          family = nb(), data = vf_final)

summary(model_gam_env_only)
AIC(model_gam_env_only)





