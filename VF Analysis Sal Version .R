# ------------------------------------------------------------
# Valley-Fever / Climate-Equity • County-Level Analysis (Final)
# ------------------------------------------------------------
rm(list = ls())

# 0. Libraries
library(tidyverse)
library(tidycensus)
library(tigris)
library(sf)
library(prism)
library(terra)
library(exactextractr)
library(lubridate)
library(here)
library(janitor)
library(stringr)
library(glmmTMB)
library(splines)

options(tigris_use_cache = TRUE, tigris_class = "sf")

# Create directories if they don't exist
if (!dir.exists(here("data_raw"))) dir.create(here("data_raw"), recursive = TRUE)
if (!dir.exists(here("data_clean"))) dir.create(here("data_clean"), recursive = TRUE)

here()
# 1. DOWNLOAD CLIMATE DATA
# ------------------------------------------------------------
# Set the download directory for PRISM data to a subfolder in data_raw
prism_set_dl_dir(here("data_raw", "prism"))

# Loop over each year from 2010 to 2023
for (year in 2010:2023) {
  tryCatch({
    # Download monthly maximum temperature (tmax) data
    get_prism_monthlys("tmax", years = year, mon = 1:12, keepZip = FALSE)
    # Download monthly precipitation (ppt) data
    get_prism_monthlys("ppt", years = year, mon = 1:12, keepZip = FALSE)
  }, error = function(e) {
    # If a year fails to download, print a message but continue
    message(paste("Could not download data for year:", year))
  })
}


# 2. EXTRACT CLIMATE DATA TO TRACTS
# ------------------------------------------------------------
# List all downloaded tmax raster files (.bil files)
pr_tmax_files <- list.files(here("data_raw", "prism"), pattern = "tmax.*\\.bil$", full.names = TRUE, recursive = TRUE)

# List all downloaded ppt raster files
pr_ppt_files  <- list.files(here("data_raw", "prism"), pattern = "ppt.*\\.bil$", full.names = TRUE, recursive = TRUE)

# Check that files exist; if not, stop execution
if (!length(pr_tmax_files) || !length(pr_ppt_files)) {
  stop("PRISM rasters missing: cannot proceed")
}

# Stack all tmax rasters into one multi-layer object
stack_tmax <- rast(pr_tmax_files)

# Stack all ppt rasters
stack_ppt  <- rast(pr_ppt_files)

# Download California census tracts shapefile (2020, cartographic boundary version)
tract_sf   <- tracts("CA", cb = TRUE, year = 2020)

# Extract mean tmax values for each tract (average over each tract area)
tmax_vals <- exact_extract(stack_tmax, tract_sf, "mean", progress = FALSE)

# Extract mean ppt values for each tract
ppt_vals  <- exact_extract(stack_ppt, tract_sf, "mean", progress = FALSE)

# Helper function to convert extracted matrix to long format with dates
munge_climate_data <- function(vals, var_name) {
  as_tibble(vals) %>%
    # Add tract ID
    mutate(tract = tract_sf$GEOID) %>%
    # Convert wide columns to long (one row per month)
    pivot_longer(cols = -tract, names_to = "original_name", values_to = var_name) %>%
    # Extract year-month from layer names and convert to date
    mutate(month = ymd(paste0(str_extract(original_name, "\\d{6}"), "01"))) %>%
    select(tract, month, {{var_name}}) %>%
    filter(!is.na(month))
}

# Create tract-level table of tmax values
climate_tract <- left_join(
  munge_climate_data(tmax_vals, "tmax"),
  munge_climate_data(ppt_vals, "ppt"),
  by = c("tract", "month")
)

# Print project directory path (for verification)
print(here())

# 3. PROCESS TRACT-LEVEL SVI & PM2.5 DATA
# ------------------------------------------------------------
# Read Social Vulnerability Index (SVI) data (CSV file in data_raw folder)
svi_path <- read_csv("data_raw/california_svi_2020.csv") 

# Read CalEnviroScreen (CES) data for California (CSV file in data_raw folder)
ces_path <- read_csv("data_raw/CalEnviroScreen_4.0_Results.csv")

# The following lines were for debugging or exploring files interactively — can be removed if no longer needed
hello <- file.choose()  # Opens file dialog to manually pick a file (not required)
glimpse(hello)          # Shows structure of that file
show(svi_path)          # Displays SVI dataframe in console
glimpse(ces_path)       # Shows structure of CES dataframe

# The file existence check below is not necessary anymore since read_csv will error if missing,
# but we keep it for reference if desired
if (!file.exists(svi_path) || !file.exists(ces_path)) stop("SVI or CES data is missing.")

# Clean and format SVI data at tract level
svi_tract <- (svi_path) %>%
  clean_names() %>%                     # Clean column names to snake_case
  filter(rpl_themes >= 0) %>%          # Keep tracts with valid overall SVI score
  transmute(
    tract       = str_pad(as.character(fips), 11, pad = "0"),  # Ensure tract codes are 11-digit strings
    svi_overall = rpl_themes                               # Overall SVI score (higher = more vulnerable)
  )

# Clean and format CES data at tract level to get PM2.5 (as dust exposure proxy)
ces_tract <- ces_path %>%
  clean_names() %>%
  transmute(
    tract = str_pad(as.character(tract), 11, pad = "0"),  # Ensure tract codes are 11-digit strings
    pm25  = pm                                           # Annual average PM2.5 concentration
  )

# The next ces_tract block is redundant and repeats the above; kept here if needed for clarity
ces_tract <- ces_path %>%
  clean_names() %>%
  transmute(
    tract = str_pad(as.character(tract), 11, pad = "0"),
    pm25  = pm
  )

# Join SVI and PM2.5 data by census tract GEOID
svi_pm_tract <- left_join(svi_tract, ces_tract, by = "tract")


# 4. PROCESS COUNTY-LEVEL VALLEY FEVER CASES
# ------------------------------------------------------------
# Read Valley Fever (VF) case data from CSV (skip first 3 header rows, provide custom column names)
vs_path <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3, 
                    col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

# Clean and standardize VF data
vf_cleaned <- vs_path %>%
  mutate(county = str_to_upper(county)) %>%    # Convert county names to uppercase for consistency
  filter(!is.na(county) &                      # Remove rows with missing county names
           !str_detect(county, "TOTAL|BERKELEY|LONG BEACH|PASADENA|\\*")) %>%  # Remove summary or city-specific rows
  mutate(
    county  = str_remove(county, " COUNTY"),   # Remove " COUNTY" suffix
    county  = str_trim(county),                # Remove any extra spaces
    cases   = cases_raw,                       # Keep raw case count
    inc_rate = parse_number(inc_rate_raw)      # Convert incidence rate string to numeric
  ) %>%
  filter(!is.na(cases) & !is.na(year) & !is.na(inc_rate)) %>%  # Keep rows with valid data
  select(county, year, cases, inc_rate)        # Keep relevant columns

# Identify "endemic" counties with ≥500 total cases and mean incidence rate >5 per 100,000
endemic_counties <- vf_cleaned %>%
  group_by(county) %>%
  summarise(
    total_cases = sum(cases, na.rm = TRUE),    # Total cases from 2001–2023
    mean_inc    = mean(inc_rate, na.rm = TRUE),# Average incidence rate
    .groups     = "drop"
  ) %>%
  filter(total_cases >= 500, mean_inc > 5) %>% # Apply thresholds
  pull(county)                                  # Return vector of qualifying county names

# Keep only endemic counties and estimate population
vf_county_cases <- vf_cleaned %>%
  filter(county %in% endemic_counties) %>% 
  mutate(population = round(cases / (inc_rate / 100000), 0)) %>%  # Back-calculate approximate population
  select(county, year, cases, population)                         # Keep final columns




# 5. BUILD COUNTY-LEVEL PANEL & RUN ANALYSIS
# ------------------------------------------------------------

# Get California county geometries from TIGER/CB shapefiles
county_data <- counties("CA", cb = TRUE, year = 2020) %>%
  as_tibble() %>%
  mutate(county = str_to_upper(NAME)) %>%    # Standardize county names to uppercase
  select(COUNTYFP, county)                   # Keep FIPS code and name

# Create tract-to-county lookup table
tract_to_county_lookup <- tracts("CA", cb = TRUE, year = 2020) %>%
  as_tibble() %>%
  left_join(county_data, by = "COUNTYFP") %>%   # Join on county FIPS code
  select(tract = GEOID, county)                 # Keep tract GEOID and county name

# Get 2020 tract-level population estimates from ACS
tract_populations <- get_acs("tract", variables = "B01003_001", state = "CA", year = 2020) %>%
  select(tract = GEOID, pop_tract = estimate)

# Aggregate tract-level SVI and PM2.5 data up to county level using population-weighted means
svi_pm_county <- svi_pm_tract %>%
  left_join(tract_to_county_lookup, by = "tract") %>%   # Add county info
  left_join(tract_populations, by = "tract") %>%       # Add tract populations for weighting
  filter(!is.na(county) & !is.na(pop_tract)) %>%
  group_by(county) %>%
  summarise(
    svi_overall = weighted.mean(svi_overall, w = pop_tract, na.rm = TRUE),  # Weighted SVI
    pm25        = weighted.mean(pm25, w = pop_tract, na.rm = TRUE)          # Weighted PM2.5
  )

# Aggregate climate data (tmax and ppt) to county-month level using mean across tracts
climate_county <- climate_tract %>%
  left_join(tract_to_county_lookup, by = "tract") %>%
  filter(!is.na(county)) %>%
  group_by(county, month) %>%
  summarise(across(c(tmax, ppt), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Expand annual county-level VF cases to monthly estimates
vf_cases_monthly <- vf_county_cases %>%
  filter(year >= 2010 & year <= 2023) %>%           # Keep years with climate data
  uncount(12, .id = "month_num") %>%               # Create 12 rows per year (1 per month)
  mutate(
    cases_monthly = round(cases / 12, 0),         # Approximate monthly cases
    month = ymd(paste(year, month_num, "01", sep = "-"))  # Create first day of each month
  ) %>%
  select(county, month, cases = cases_monthly, population)

# Join VF data with county-level SVI/PM2.5 and climate data to create full panel
county_month_panel <- vf_cases_monthly %>%
  left_join(svi_pm_county, by = "county") %>%
  left_join(climate_county, by = c("county", "month")) %>%
  drop_na()  # Drop rows with missing data

# --- Statistical Model (Publication Version) ---
# Prepare data with time variable and factors
model_data <- county_month_panel %>%
  mutate(
    time = as.numeric(month - min(month)) / 30.25,   # Approximate months since start (scaled)
    month_factor = factor(month(month)),             # Month categorical variable for seasonality
    county = factor(county)                          # County as a factor
  )

# Final publication-ready model including PM2.5 as a dust control
final_model_pub <- glmmTMB(
  cases ~ svi_overall + pm25 + tmax + ppt +     # Predictors: SVI, PM2.5, temp, precipitation
    ns(time, df = 4) +                          # Smooth term for time trend (4 degrees freedom)
    month_factor +                              # Seasonal (monthly) effect
    factor(county) +                            # County fixed effects
    offset(log(population)),                    # Offset for log(population)
  family = nbinom2,                             # Negative binomial model
  data = model_data
)

print("--- Final Model Results ---")
print(summary(final_model_pub))  # Show model summary with coefficients

# --- Final Result Interpretation ---
print("--- SVI Incidence Rate Ratio (IRR) with 95% CI ---")

# Get confidence intervals for coefficients
model_ci <- confint(final_model_pub, parm = "beta_", method = "wald")
model_coeffs <- summary(final_model_pub)$coefficients$cond

# Extract log(IRR) and CI for SVI
svi_log_irr <- model_coeffs["svi_overall", "Estimate"]
svi_ci_log <- model_ci["svi_overall", ]

# Compute IRR and CI for a 0.1 unit increase in SVI
irr_0.1_svi <- exp(svi_log_irr * 0.1)
ci_0.1_svi <- exp(svi_ci_log * 0.1)

cat(sprintf(
  "IRR for a 0.1 unit increase in SVI: %.3f (95%% CI: %.3f, %.3f)\n",
  irr_0.1_svi,
  ci_0.1_svi[1],
  ci_0.1_svi[2]
))


# -------------------------------------------------------------------
# ADD-ON: Re-run the model WITHOUT the dust confounder for comparison
# -------------------------------------------------------------------
print("--- Running Model WITHOUT Dust Control (for comparison) ---")

# Define model without PM2.5 term to assess sensitivity
model_without_dust <- glmmTMB(
  cases ~ svi_overall + tmax + ppt +
    ns(time, df = 4) +
    month_factor +
    factor(county) +
    offset(log(population)),
  family = nbinom2,
  data = model_data
)

print(summary(model_without_dust))  # Show summary

# --- Interpretation for the model without dust control ---
print("--- SVI IRR (95% CI) WITHOUT Dust Control ---")

# Extract CI and coefficients
model_ci_no_dust <- confint(model_without_dust, parm = "beta_", method = "wald")
model_coeffs_no_dust <- summary(model_without_dust)$coefficients$cond

svi_log_irr_no_dust <- model_coeffs_no_dust["svi_overall", "Estimate"]
svi_ci_log_no_dust <- model_ci_no_dust["svi_overall", ]

# Compute IRR and CI for a 0.1 unit increase in SVI
irr_0.1_svi_no_dust <- exp(svi_log_irr_no_dust * 0.1)
ci_0.1_svi_no_dust <- exp(svi_ci_log_no_dust * 0.1)

cat(sprintf(
  "IRR for a 0.1 unit increase in SVI (no dust control): %.3f (95%% CI: %.3f, %.3f)\n",
  irr_0.1_svi_no_dust,
  ci_0.1_svi_no_dust[1],
  ci_0.1_svi_no_dust[2]
))


# -------------------------------------------------------------------
# SENSITIVITY TEST: EXCLUDING URBAN TRACTS
# -------------------------------------------------------------------
print("--- Running Sensitivity Test: Excluding Urban Tracts ---")

# 1. Identify non-urban tracts based on paper's criteria
#    (Keep tracts BELOW the 40th percentile of population density)

# Get tract land areas from TIGER/CB shapefile
tract_areas <- tracts("CA", cb = TRUE, year = 2020) %>%
  as_tibble() %>%
  select(tract = GEOID, area_sq_m = ALAND)  # ALAND = land area in square meters

# Join tract population estimates to areas to calculate density
pop_density <- left_join(tract_populations, tract_areas, by = "tract") %>%
  mutate(
    density = pop_tract / (area_sq_m / 1e6)  # People per square kilometer (convert m² to km²)
  ) %>%
  filter(!is.na(density) & is.finite(density))  # Keep valid density rows

# Find 40th percentile density value (threshold for defining "urban")
urban_threshold <- quantile(pop_density$density, 0.40, na.rm = TRUE)

# Create list of non-urban tracts (tracts at or below the threshold)
non_urban_tracts <- pop_density %>%
  filter(density <= urban_threshold) %>%
  pull(tract)

# 2. Re-aggregate SVI and PM2.5 using only non-urban tracts
svi_pm_county_no_urban <- svi_pm_tract %>%
  filter(tract %in% non_urban_tracts) %>%         # Keep only non-urban tracts
  left_join(tract_to_county_lookup, by = "tract") %>%
  left_join(tract_populations, by = "tract") %>%
  filter(!is.na(county) & !is.na(pop_tract)) %>%
  group_by(county) %>%
  summarise(
    svi_overall = weighted.mean(svi_overall, w = pop_tract, na.rm = TRUE),  # Population-weighted SVI
    pm25        = weighted.mean(pm25, w = pop_tract, na.rm = TRUE)          # Population-weighted PM2.5
  )

# Aggregate climate data to county level using only non-urban tracts
climate_county_no_urban <- climate_tract %>%
  filter(tract %in% non_urban_tracts) %>%         # Keep only non-urban tracts
  left_join(tract_to_county_lookup, by = "tract") %>%
  filter(!is.na(county)) %>%
  group_by(county, month) %>%
  summarise(across(c(tmax, ppt), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# 3. Build new county-month panel using filtered data
county_month_panel_no_urban <- vf_cases_monthly %>%
  left_join(svi_pm_county_no_urban, by = "county") %>%
  left_join(climate_county_no_urban, by = c("county", "month")) %>%
  drop_na()  # Remove rows with missing data

# 4. Re-run the model using data excluding urban tracts
model_data_no_urban <- county_month_panel_no_urban %>%
  mutate(
    time = as.numeric(month - min(month)) / 30.25,  # Approximate months since start
    month_factor = factor(month(month)),            # Month factor for seasonality
    county = factor(county)                         # County factor
  )

# Define and fit the negative binomial model
final_model_sens_urban <- glmmTMB(
  cases ~ svi_overall + pm25 + tmax + ppt + 
    ns(time, df = 4) + 
    month_factor + 
    factor(county) + 
    offset(log(population)),
  family = nbinom2,
  data = model_data_no_urban
)

# Print results for the sensitivity model
print("--- Sensitivity Model Results (Excluding Urban Tracts) ---")
print(summary(final_model_sens_urban))

# Install DHARMa package if needed for residual checks (diagnostics)
install.packages("DHARMa")

# -------------------------------------------------------------------
# MODEL DIAGNOSTICS WITH DHARMa
# -------------------------------------------------------------------
print("--- Running DHARMa Model Diagnostics ---")

# Load DHARMa package for residual diagnostics
# Only run install.packages("DHARMa") if you haven't already installed it
# install.packages("DHARMa")
library(DHARMa)
library(splines)  # Already loaded earlier, repeated here for clarity

# Simulate residuals from your final model (with dust control)
# n = number of simulations; 250 is a standard choice
simulation_output <- simulateResiduals(fittedModel = final_model_pub, n = 250)

# Create standard DHARMa diagnostic plots
# This will open a plot window (if using RStudio or similar)
plot(simulation_output)

# Test for overdispersion (checks if variance matches model assumptions)
testDispersion(simulation_output)

# Test for zero-inflation (checks for excess zeros not captured by model)
zi <- testZeroInflation(simulation_output)
print(zi)


# -------------------------------------------------------------------
# MODEL COMPARISON: ZINB and NB1
# -------------------------------------------------------------------
library(glmmTMB)
print("--- Running Additional Models for Comparison ---")

# 1. Zero-Inflated Negative Binomial (nbinom2) model
# Adds a zero-inflation component to account for extra zeros
print("--- Fitting Zero-Inflated NB2 Model (ZINB) ---")
model_zinb <- glmmTMB(
  cases ~ svi_overall + pm25 + tmax + ppt +
    ns(time, df = 4) +
    month_factor +
    factor(county) +
    offset(log(population)),
  ziformula = ~1,      # Add intercept-only zero-inflation component
  family = nbinom2,    # Negative binomial (variance = mu + mu^2 / theta)
  data = model_data
)
print(summary(model_zinb))

# 2. Standard Negative Binomial 1 (nbinom1) model
# Alternative parameterization: variance = phi * mu
print("--- Fitting NB1 Model ---")
model_nb1 <- glmmTMB(
  cases ~ svi_overall + pm25 + tmax + ppt +
    ns(time, df = 4) +
    month_factor +
    factor(county) +
    offset(log(population)),
  family = nbinom1,  # Different variance structure
  data = model_data
)
print(summary(model_nb1))

# 3. Compare models using AIC
# Lower AIC indicates a better trade-off between fit and complexity
print("--- Model Comparison (AIC) ---")

# Make sure your original negative binomial model is assigned correctly
# Here we assume 'final_model_fe' is your standard NB2 model
# Adjust if your object name is different
final_model_nb2 <- final_model_pub


print(AIC(final_model_nb2, model_zinb, model_nb1))


# -------------------------------------------------------------------
# CHECKING FOR SPATIAL AUTOCORRELATION (Corrected for Repeated Measures)
# -------------------------------------------------------------------
print("--- Running Spatial Autocorrelation Test ---")

# 1. Simulate residuals again (here from 'final_model_fe')
simulation_output <- simulateResiduals(fittedModel = final_model_pub, n = 250)

# 2. Aggregate residuals by county
# This accounts for repeated measurements by averaging within counties
aggregated_residuals <- recalculateResiduals(
  simulationOutput = simulation_output,
  group = model_data$county
)

# 3. Get unique county centroids (approximate coordinates)
unique_county_coords <- model_data %>%
  distinct(county) %>%
  left_join(
    counties("CA", cb = TRUE, year = 2020) %>%
      mutate(county = str_to_upper(NAME)) %>%
      st_centroid() %>%                                           # Calculate county centroids
      mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
      as_tibble() %>%
      select(county, lon, lat),
    by = "county"
  )

# 4. Run DHARMa spatial autocorrelation test (Moran's I test)
spatial_test_result <- testSpatialAutocorrelation(
  simulationOutput = aggregated_residuals,
  x = unique_county_coords$lon,
  y = unique_county_coords$lat,
  plot = TRUE   # Create a diagnostic plot
)

print("--- Moran's I Test for Spatial Autocorrelation ---")
print(spatial_test_result)


################################### Graph #####################################

# --- Fit simple linear regression ---
model_linear_simple <- lm(cases ~ svi_overall, data = model_data)
summary_model <- summary(model_linear_simple)

# --- Extract slope coefficient and p-value ---
coef_svi <- summary_model$coefficients["svi_overall", "Estimate"]
pval_svi <- summary_model$coefficients["svi_overall", "Pr(>|t|)"]

# --- Format annotation text ---
annotation_text <- sprintf("Slope: %.1f\np-value: %.3f", coef_svi, pval_svi)

# --- Create plot with annotation ---
library(ggplot2)

ggplot(model_data, aes(x = svi_overall, y = cases)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Valley Fever Cases vs. SVI (County-Month Panel, Simple Linear Regression)",
    x = "SVI (Overall)",
    y = "Monthly VF Cases"
  ) +
  annotate("text", x = min(model_data$svi_overall, na.rm = TRUE) + 0.02, 
           y = max(model_data$cases, na.rm = TRUE) * 0.9,
           label = annotation_text, hjust = 0, size = 4, color = "black") +
  theme_minimal()

install.packages("broom")
library(broom)

tidy(model_linear_basic, conf.int = TRUE)


library(tidyverse)

# Read Valley Fever data
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

# Read SVI data
svi <- read_csv("data_raw/california_svi_2020.csv")

vf_clean <- vf_cases %>%
  mutate(county = str_to_upper(county)) %>%
  filter(!is.na(county) & !str_detect(county, "TOTAL|BERKELEY|LONG BEACH|PASADENA|\\*")) %>%
  mutate(county = str_remove(county, " COUNTY"),
         county = str_trim(county),
         cases = as.numeric(cases_raw)) %>%
  filter(!is.na(cases)) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop")

svi_clean <- svi %>%
  janitor::clean_names() %>%
  filter(rpl_themes >= 0) %>%
  transmute(county = str_to_upper(county),
            svi_overall = rpl_themes)

vf_svi <- left_join(vf_clean, svi_clean, by = "county") %>%
  drop_na()

ggplot(vf_svi, aes(x = svi_overall, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Total Valley Fever Cases vs. SVI (County-level)",
    x = "SVI (Overall)",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()


# VF counties
vf_clean %>% distinct(county) %>% count()

# SVI counties
svi_clean %>% distinct(county) %>% count()

ggplot(vf_svi, aes(x = svi_overall, y = total_cases)) +
  geom_point(size = 2, alpha = 0.6, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Total Valley Fever Cases vs. SVI (County-level)",
    x = "SVI (Overall)",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()

svi %>% janitor::clean_names() %>% distinct(county) %>% count()

svi %>% janitor::clean_names() %>% count(county)

svi_clean <- svi %>%
  janitor::clean_names() %>%
  filter(rpl_themes >= 0) %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(svi_overall = mean(rpl_themes, na.rm = TRUE), .groups = "drop")

vf_svi <- left_join(vf_clean, svi_clean, by = "county") %>%
  drop_na()

nrow(vf_svi)  # Should now be ~58

vf_svi %>% count(county)

ggplot(vf_svi, aes(x = svi_overall, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Total Valley Fever Cases vs. SVI (County-level)",
    x = "SVI (Overall, aggregated by county)",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()

counties_over_500 <- vf_svi %>%
  filter(total_cases > 500) %>%
  pull(county)

vf_svi_filtered <- vf_svi %>%
  filter(county %in% counties_over_500)

model_lm_filtered <- lm(total_cases ~ svi_overall, data = vf_svi_filtered)
summary(model_lm_filtered)

library(ggplot2)

ggplot(vf_svi_filtered, aes(x = svi_overall, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Total VF Cases vs. SVI (Counties >500 Cases)",
    x = "SVI (Overall)",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()


vf_cases_clean <- vf_cases %>%
  mutate(county = str_to_upper(county)) %>%
  mutate(county = str_remove(county, " COUNTY"),
         county = str_trim(county),
         cases = as.numeric(cases_raw)) %>%
  filter(!is.na(cases))

special_cities <- vf_cases_clean %>%
  filter(county %in% c("BERKELEY", "LONG BEACH", "PASADENA")) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop")

vf_cases_summary <- vf_cases_clean %>%
  mutate(county = case_when(
    county == "BERKELEY" ~ "ALAMEDA",
    county %in% c("LONG BEACH", "PASADENA") ~ "LOS ANGELES",
    TRUE ~ county
  )) %>%
  group_by(county, year) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
  group_by(county) %>%
  summarise(
    total_cases = sum(cases, na.rm = TRUE),
    avg_cases_per_year = mean(cases, na.rm = TRUE),
    .groups = "drop"
  )

eligible_counties <- vf_cases_summary %>%
  filter(total_cases > 500, avg_cases_per_year >= 5) %>%
  pull(county)

vf_svi_filtered <- vf_svi %>%
  filter(county %in% eligible_counties)

model_lm_filtered <- lm(total_cases ~ svi_overall, data = vf_svi_filtered)
summary(model_lm_filtered)

ggplot(vf_svi_filtered, aes(x = svi_overall, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Total VF Cases vs. SVI (Corrected Counties >500 & ≥5/year)",
    x = "SVI (Overall)",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()

svi <- read_csv("data_raw/california_svi_2020.csv") %>%
  clean_names()

names(svi)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
library(janitor)

# Read & clean
svi <- read_csv("data_raw/california_svi_2020.csv")
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(county = str_to_upper(county)) %>%
  # Remove only TOTAL and "*" rows
  filter(!is.na(county) & !str_detect(county, "TOTAL|\\*")) %>%
  mutate(county = str_remove(county, " COUNTY"),
         county = str_trim(county),
         cases = as.numeric(cases_raw),
         # Reassign
         county = case_when(
           county == "BERKELEY" ~ "ALAMEDA",
           county == "LONG BEACH" ~ "LOS ANGELES",
           county == "PASADENA" ~ "LOS ANGELES",
           TRUE ~ county
         )) %>%
  filter(!is.na(cases)) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop")

svi_clean <- svi %>%
  janitor::clean_names() %>%
  filter(rpl_themes >= 0) %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    pov = mean(ep_pov150, na.rm = TRUE),
    unemp = mean(ep_unemp, na.rm = TRUE),
    hburd = mean(ep_hburd, na.rm = TRUE),
    nohsdp = mean(ep_nohsdp, na.rm = TRUE),
    uninsur = mean(ep_uninsur, na.rm = TRUE),
    .groups = "drop"
  )

vf_svi <- left_join(vf_clean, svi_clean, by = "county") %>%
  drop_na()

# Pivot longer
vf_long <- vf_svi %>%
  pivot_longer(cols = c(pov, unemp, hburd, nohsdp, uninsur), 
               names_to = "variable", 
               values_to = "svi_value")

# Facet plot
ggplot(vf_long, aes(x = svi_value, y = total_cases)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Valley Fever Cases vs. Theme 1 SVI Variables (County-level)",
    x = "SVI Variable Percentile Score",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()

# Individual regressions
theme1_vars <- c("pov", "unemp", "hburd", "nohsdp", "uninsur")

for (var in theme1_vars) {
  df <- vf_svi %>%
    select(county, total_cases, !!sym(var)) %>%
    filter(!is.na(total_cases), !is.na(!!sym(var))) %>%
    rename(var_value = !!sym(var))
  
  model <- lm(total_cases ~ var_value, data = df)
  
  cat("\n\n----------------------------------------\n")
  cat("Summary for variable:", var, "\n")
  cat("----------------------------------------\n")
  print(summary(model))
}

#################### Weakest predictors unemp, uninsur #######################




# Multiple regression with all Theme 1 variables
df_theme1 <- vf_svi %>%
  select(county, total_cases, pov, unemp, hburd, nohsdp, uninsur) %>%
  filter(
    !is.na(total_cases), 
    !is.na(pov), 
    !is.na(unemp), 
    !is.na(hburd), 
    !is.na(nohsdp), 
    !is.na(uninsur)
  )

model_theme1 <- lm(total_cases ~ pov + unemp + hburd + nohsdp + uninsur, data = df_theme1)

# Print full summary
summary(model_theme1)




# Read & clean
svi <- read_csv("data_raw/california_svi_2020.csv")
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(county = str_to_upper(county)) %>%
  filter(!is.na(county) & !str_detect(county, "TOTAL|\\*")) %>%
  mutate(county = str_remove(county, " COUNTY"),
         county = str_trim(county),
         cases = as.numeric(cases_raw),
         county = case_when(
           county == "BERKELEY" ~ "ALAMEDA",
           county == "LONG BEACH" ~ "LOS ANGELES",
           county == "PASADENA" ~ "LOS ANGELES",
           TRUE ~ county
         )) %>%
  filter(!is.na(cases)) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop")

svi_clean <- svi %>%
  janitor::clean_names() %>%
  filter(rpl_themes >= 0) %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    age65 = mean(ep_age65, na.rm = TRUE),
    age17 = mean(ep_age17, na.rm = TRUE),
    disabl = mean(ep_disabl, na.rm = TRUE),
    sngpnt = mean(ep_sngpnt, na.rm = TRUE),
    .groups = "drop"
  )

vf_svi <- left_join(vf_clean, svi_clean, by = "county") %>%
  drop_na()

vf_long <- vf_svi %>%
  pivot_longer(cols = c(age65, age17, disabl, sngpnt), 
               names_to = "variable", 
               values_to = "svi_value")

ggplot(vf_long, aes(x = svi_value, y = total_cases)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Valley Fever Cases vs. Theme 2 SVI Variables (County-level)",
    x = "SVI Variable Percentile Score",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()

theme2_vars <- c("age65", "age17", "disabl", "sngpnt")

for (var in theme2_vars) {
  df <- vf_svi %>%
    select(county, total_cases, !!sym(var)) %>%
    filter(!is.na(total_cases), !is.na(!!sym(var))) %>%
    rename(var_value = !!sym(var))
  
  model <- lm(total_cases ~ var_value, data = df)
  
  cat("\n\n----------------------------------------\n")
  cat("Summary for variable:", var, "\n")
  cat("----------------------------------------\n")
  print(summary(model))
}

################# Weakest predictor: disabl and neg for age65 ##################

df_theme2 <- vf_svi %>%
  select(county, total_cases, age65, age17, disabl, sngpnt) %>%
  filter(
    !is.na(total_cases), 
    !is.na(age65), 
    !is.na(age17), 
    !is.na(disabl), 
    !is.na(sngpnt)
  )

model_theme2 <- lm(total_cases ~ age65 + age17 + disabl + sngpnt, data = df_theme2)

summary(model_theme2) #disabl is the only negative





svi_clean <- svi %>%
  janitor::clean_names() %>%
  filter(rpl_themes >= 0) %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    minrty = mean(ep_minrty, na.rm = TRUE),
    limeng = mean(ep_limeng, na.rm = TRUE),
    .groups = "drop"
  )

vf_svi <- left_join(vf_clean, svi_clean, by = "county") %>%
  drop_na()

vf_long <- vf_svi %>%
  pivot_longer(cols = c(minrty, limeng), 
               names_to = "variable", 
               values_to = "svi_value")

ggplot(vf_long, aes(x = svi_value, y = total_cases)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Valley Fever Cases vs. Theme 3 SVI Variables (County-level)",
    x = "SVI Variable Percentile Score",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()

theme3_vars <- c("minrty", "limeng")

for (var in theme3_vars) {
  df <- vf_svi %>%
    select(county, total_cases, !!sym(var)) %>%
    filter(!is.na(total_cases), !is.na(!!sym(var))) %>%
    rename(var_value = !!sym(var))
  
  model <- lm(total_cases ~ var_value, data = df)
  
  cat("\n\n----------------------------------------\n")
  cat("Summary for variable:", var, "\n")
  cat("----------------------------------------\n")
  print(summary(model))
}

df_theme3 <- vf_svi %>%
  select(county, total_cases, minrty, limeng) %>%
  filter(
    !is.na(total_cases), 
    !is.na(minrty), 
    !is.na(limeng)
  )

model_theme3 <- lm(total_cases ~ minrty + limeng, data = df_theme3)

summary(model_theme3) #limeng has nothing 




svi_clean <- svi %>%
  janitor::clean_names() %>%
  filter(rpl_themes >= 0) %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    munit = mean(ep_munit, na.rm = TRUE),
    mobile = mean(ep_mobile, na.rm = TRUE),
    crowd = mean(ep_crowd, na.rm = TRUE),
    noveh = mean(ep_noveh, na.rm = TRUE),
    groupq = mean(ep_groupq, na.rm = TRUE),
    .groups = "drop"
  )

vf_svi <- left_join(vf_clean, svi_clean, by = "county") %>%
  drop_na()

vf_long <- vf_svi %>%
  pivot_longer(cols = c(munit, mobile, crowd, noveh, groupq), 
               names_to = "variable", 
               values_to = "svi_value")

ggplot(vf_long, aes(x = svi_value, y = total_cases)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Valley Fever Cases vs. Theme 4 SVI Variables (County-level)",
    x = "SVI Variable Percentile Score",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()

theme4_vars <- c("munit", "mobile", "crowd", "noveh", "groupq")

for (var in theme4_vars) {
  df <- vf_svi %>%
    select(county, total_cases, !!sym(var)) %>%
    filter(!is.na(total_cases), !is.na(!!sym(var))) %>%
    rename(var_value = !!sym(var))
  
  model <- lm(total_cases ~ var_value, data = df)
  
  cat("\n\n----------------------------------------\n")
  cat("Summary for variable:", var, "\n")
  cat("----------------------------------------\n")
  print(summary(model))
}

######################munit, mobile, noveh, groupq — all non-significant.#####

df_theme4 <- vf_svi %>%
  select(county, total_cases, munit, mobile, crowd, noveh, groupq) %>%
  filter(
    !is.na(total_cases), 
    !is.na(munit), 
    !is.na(mobile), 
    !is.na(crowd), 
    !is.na(noveh), 
    !is.na(groupq)
  )

model_theme4 <- lm(total_cases ~ munit + mobile + crowd + noveh + groupq, data = df_theme4)

summary(model_theme4)


library(tidyverse)
library(janitor)

# --- Read data ---
svi <- read_csv("data_raw/california_svi_2020.csv")
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

# --- Clean Valley Fever data ---
vf_clean <- vf_cases %>%
  mutate(county = str_to_upper(county)) %>%
  filter(!is.na(county) & !str_detect(county, "TOTAL|\\*")) %>%
  mutate(county = str_remove(county, " COUNTY"),
         county = str_trim(county),
         cases = as.numeric(cases_raw),
         county = case_when(
           county == "BERKELEY" ~ "ALAMEDA",
           county == "LONG BEACH" ~ "LOS ANGELES",
           county == "PASADENA" ~ "LOS ANGELES",
           TRUE ~ county
         )) %>%
  filter(!is.na(cases)) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop")

# --- Clean SVI data ---
svi_clean <- svi %>%
  janitor::clean_names() %>%
  filter(rpl_themes >= 0) %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(svi_overall = mean(rpl_themes, na.rm = TRUE), .groups = "drop")

# --- Merge ---
vf_svi <- left_join(vf_clean, svi_clean, by = "county") %>%
  drop_na()

# --- Plot ---
ggplot(vf_svi, aes(x = svi_overall, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Total Valley Fever Cases vs. Overall SVI (All Counties)",
    x = "SVI Overall Score",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()

# --- Model & summary ---
model_svi <- lm(total_cases ~ svi_overall, data = vf_svi)
summary(model_svi)



library(tidyverse)
library(janitor)

library(tidyverse)

# --- Read Valley Fever data ---
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(county = str_to_upper(county)) %>%
  filter(!is.na(county) & !str_detect(county, "TOTAL|\\*")) %>%
  mutate(county = str_remove(county, " COUNTY"),
         county = str_trim(county),
         cases = as.numeric(cases_raw),
         county = case_when(
           county == "BERKELEY" ~ "ALAMEDA",
           county == "LONG BEACH" ~ "LOS ANGELES",
           county == "PASADENA" ~ "LOS ANGELES",
           TRUE ~ county
         )) %>%
  filter(!is.na(cases)) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop")

# --- Read SVI data ---
svi <- read_csv("data_raw/california_svi_2020.csv")

# --- Clean SVI and compute custom raw score ---
svi_clean <- svi %>%
  clean_names() %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    pov = mean(ep_pov150, na.rm = TRUE),
    hburd = mean(ep_hburd, na.rm = TRUE),
    nohsdp = mean(ep_nohsdp, na.rm = TRUE),
    age65 = mean(ep_age65, na.rm = TRUE),
    age17 = mean(ep_age17, na.rm = TRUE),
    sngpnt = mean(ep_sngpnt, na.rm = TRUE),
    minrty = mean(ep_minrty, na.rm = TRUE),
    crowd = mean(ep_crowd, na.rm = TRUE),
    limeng = mean(ep_limeng, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(custom_svi_raw = mean(c(pov, hburd, nohsdp, age65, age17, sngpnt, minrty, crowd, limeng), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(custom_svi_percentile = percent_rank(custom_svi_raw)) # ✅ HERE we convert to 0-1 percentile

# --- Merge VF with custom SVI ---
vf_svi <- left_join(vf_clean, svi_clean, by = "county") %>%
  drop_na()

# --- Plot custom percentile SVI vs total VF cases ---
ggplot(vf_svi, aes(x = custom_svi_percentile, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Total Valley Fever Cases vs. Custom SVI (All Counties, Percentile Scale)",
    x = "Custom SVI Percentile Score (0–1)",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()

# --- Statistical summary ---
model_custom_svi <- lm(total_cases ~ custom_svi_percentile, data = vf_svi)
summary(model_custom_svi)






library(MASS)

# --- Read Valley Fever data ---
vf_cases <- read_csv("valley_fever_cases_by_lhd_2001-2023.csv", skip = 3,
                     col_names = c("county", "year", "cases_raw", "inc_rate_raw"))

vf_clean <- vf_cases %>%
  mutate(county = str_to_upper(county)) %>%
  filter(!is.na(county) & !str_detect(county, "TOTAL|\\*")) %>%
  mutate(county = str_remove(county, " COUNTY"),
         county = str_trim(county),
         cases = as.numeric(cases_raw),
         county = case_when(
           county == "BERKELEY" ~ "ALAMEDA",
           county == "LONG BEACH" ~ "LOS ANGELES",
           county == "PASADENA" ~ "LOS ANGELES",
           TRUE ~ county
         )) %>%
  filter(!is.na(cases)) %>%
  group_by(county) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE), .groups = "drop")

# --- Read SVI data ---
svi <- read_csv("data_raw/california_svi_2020.csv")

# --- Clean SVI and compute custom SVI ---
svi_clean <- svi %>%
  clean_names() %>%
  group_by(county = str_to_upper(county)) %>%
  summarise(
    pov = mean(ep_pov150, na.rm = TRUE),
    hburd = mean(ep_hburd, na.rm = TRUE),
    nohsdp = mean(ep_nohsdp, na.rm = TRUE),
    age65 = mean(ep_age65, na.rm = TRUE),
    age17 = mean(ep_age17, na.rm = TRUE),
    sngpnt = mean(ep_sngpnt, na.rm = TRUE),
    minrty = mean(ep_minrty, na.rm = TRUE),
    crowd = mean(ep_crowd, na.rm = TRUE),
    limeng = mean(ep_limeng, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(custom_svi = mean(c(pov, hburd, nohsdp, age65, age17, sngpnt, minrty, crowd, limeng), na.rm = TRUE)) %>%
  ungroup()

# --- Scale custom SVI to 0–1 percentile scale ---
svi_clean <- svi_clean %>%
  mutate(custom_svi_percentile = percent_rank(custom_svi))

# --- Merge VF with custom SVI ---
vf_svi <- left_join(vf_clean, svi_clean, by = "county") %>%
  drop_na()

# --- Fit negative binomial model using percentile scaled SVI ---
nb_model <- glm.nb(total_cases ~ custom_svi_percentile, data = vf_svi)

# Create new data frame for predictions
new_data <- data.frame(custom_svi_percentile = seq(0, 1, length.out = 100))
new_data$pred_cases <- predict(nb_model, newdata = new_data, type = "response")

# --- Plot ---
ggplot(vf_svi, aes(x = custom_svi_percentile, y = total_cases)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_line(data = new_data, aes(x = custom_svi_percentile, y = pred_cases), color = "blue", size = 1.2) +
  labs(
    title = "Negative Binomial Model: Valley Fever Cases vs. Custom SVI (0–1 Scale)",
    x = "Custom SVI Percentile (0–1)",
    y = "Total VF Cases (2001–2023)"
  ) +
  theme_minimal()

# --- Print summary ---
summary(nb_model)
