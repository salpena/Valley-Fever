# Clear all objects from the workspace
rm(list = ls())

# Install packages (only once)
install.packages("tidycensus")
install.packages("zctaCrosswalk")
install.packages("readxl")
install.packages("tigris")
install.packages("sf")
install.packages("leaflet")
install.packages("psych")
install.packages("mgcv")
install.packages("viridis")

# Load all libraries
library(tidyverse)    # Collection of data tools
library(ggplot2)      # For plotting
library(readr)        # CSV reading
library(zctaCrosswalk)# ZCTA tools
library(tidycensus)   # Census API
library(dplyr)        # Data manipulation
library(readxl)       # Excel reading
library(tigris)       # Shapefiles
library(sf)           # Spatial data
library(leaflet)      # Interactive maps
library(psych)        # Phi coefficient
library(mgcv)         # GAM modeling
library(viridis)      # Color scales
library(stringr)      # String manipulation



# Check tidycensus package (just shows its content)
head(tidycensus)

# Set Census API key (replace with your own key if needed)
census_api_key("94dc655738967530ffcdc4dc68b34b420b902d6a", install = TRUE)

# Reload environment so the key is active immediately
readRenviron("~/.Renviron")  

# Define the variables you want to pull from ACS
vars <- c(
  income = "B19013_001",       # Median household income
  pop_total = "B01003_001",    # Total population
  white = "B02001_002",        # White alone
  black = "B02001_003",        # Black or African American alone
  hispanic = "B03003_003",     # Hispanic or Latino
  poverty = "B17001_002"       # Population below poverty line
)

# Get ACS data at the ZCTA level for 2022
zcta_data <- get_acs(
  geography = "zcta",
  variables = vars,
  year = 2022,
  survey = "acs5",
  output = "wide"
)

# Get list of California ZCTAs
ca_zctas <- get_zctas_by_state("CA")

# Filter ACS data to only include California ZCTAs
ca_data <- zcta_data %>%
  filter(GEOID %in% ca_zctas)

# Select and rename relevant columns
ca_data_clean <- ca_data %>%
  select(
    zcta = GEOID,
    income = incomeE,
    total_pop = pop_totalE,
    white = whiteE,
    black = blackE,
    hispanic = hispanicE,
    poverty = povertyE
  )

# Create new percentage columns for demographics and poverty
ca_data_clean <- ca_data_clean %>%
  mutate(
    pct_white = white / total_pop * 100,
    pct_black = black / total_pop * 100,
    pct_hispanic = hispanic / total_pop * 100,
    pct_poverty = poverty / total_pop * 100
  )

# Print first 10 rows in console
head(ca_data_clean, 10)

# Or open a spreadsheet-like view in RStudio
View(ca_data_clean)

# Read CalEnviroScreen data from CSV
ces_data <- read_csv("CalEnviroScreen_4.0_Results.csv")
View(ces_data)  # View entire dataset

# Check column names in CES data
colnames(ces_data)

# Select and rename relevant columns from CES data
ces_clean <- ces_data %>%
  select(
    tract = tract,
    total_pop = ACS2019TotalPop,
    pm25 = pm,
    pm25_percentile = pmP,
    ces_score = CIscore,
    ces_score_percentile = CIscoreP,
    pct_hispanic = Hispanic_pct,
    pct_white = White_pct,
    pct_black = African_American_pct,
    pct_poverty = pov
  )

# Check first rows
head(ces_clean)

# Read ZIP-to-tract crosswalk Excel file
crosswalk <- read_excel("ZIP_TRACT_122022.xlsx")

# Check column names and first rows
colnames(crosswalk)
head(crosswalk)

# Filter crosswalk to only California and prepare columns
crosswalk_ca <- crosswalk %>%
  filter(USPS_ZIP_PREF_STATE == "CA") %>%
  mutate(
    tract = as.character(TRACT),
    zcta = as.character(ZIP)
  ) %>%
  select(tract, zcta)

# Check tract values in crosswalk and CES data
head(crosswalk_ca$tract)
head(ces_clean$tract) #numbers should not have quotations 

# Format tract IDs to have 11 characters
crosswalk_ca <- crosswalk_ca %>%
  mutate(tract = sprintf("%011s", tract))

ces_clean <- ces_clean %>%
  mutate(tract = sprintf("%011s", tract))

# Join crosswalk and CES data by tract
ces_zip <- left_join(crosswalk_ca, ces_clean, by = "tract")

# Summary statistics of PM2.5
summary(ces_zip$pm25)

# Check joined data
head(ces_zip)

# Group CES data by ZIP and calculate average metrics
ces_zip_summary <- ces_zip %>%
  group_by(zcta) %>%
  summarize(
    avg_pm25 = mean(pm25, na.rm = TRUE),
    avg_ces_score = mean(ces_score, na.rm = TRUE),
    avg_pct_poverty = mean(pct_poverty, na.rm = TRUE),
    avg_pct_hispanic = mean(pct_hispanic, na.rm = TRUE),
    .groups = "drop"
  )

# Join ACS data with CES summary data by ZIP
ces_full <- inner_join(ca_data_clean, ces_zip_summary, by = "zcta")

# Create bar plot of top 20 ZIP codes by PM2.5
ces_full %>%
  arrange(desc(avg_pm25)) %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = reorder(zcta, avg_pm25), y = avg_pm25)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(title = "Top 20 Most Polluted ZIP Codes (PM2.5)",
       x = "ZIP Code", y = "Average PM2.5")

# Enable caching for shapefiles
options(tigris_use_cache = TRUE)

# Download ZCTA shapefiles
zcta_shapes <- zctas(cb = TRUE, year = 2020)
colnames(zcta_shapes)

# Optionally filter for ZIP codes starting with "9" (California)
zcta_shapes <- zctas(cb = TRUE, year = 2020) %>%
  filter(startsWith(ZCTA5CE20, "9"))

# Filter shapes to ZIP codes present in CES data
zcta_shapes_ca <- zcta_shapes %>%
  filter(ZCTA5CE20 %in% ces_zip_summary$zcta)

print(missing_shapes$zcta)


# Convert ZIP codes to character for merging
ces_zip_summary <- ces_zip_summary %>%
  mutate(zcta = as.character(zcta))

zcta_shapes_ca <- zcta_shapes_ca %>%
  mutate(zcta = as.character(ZCTA5CE20))

# Join shapefile data with CES summary data
ces_map_data <- left_join(zcta_shapes_ca, ces_zip_summary, by = "zcta")

# Plot static map with ggplot
ggplot(data = ces_map_data) +
  geom_sf(aes(fill = avg_pm25), color = NA) +
  scale_fill_viridis_c(option = "C", na.value = "grey90", name = "PM2.5") +
  labs(
    title = "Air Pollution (PM2.5) by ZIP Code in California",
    subtitle = "Using CalEnviroScreen 4.0 + ZIP Crosswalk",
    caption = "ZIPs with missing data shown in grey"
  ) +
  theme_minimal()

# Create interactive map
leaflet(data = ces_map_data) %>%
  addProviderTiles("CartoDB.Positron") %>%  # Base map
  addPolygons(
    fillColor = ~colorQuantile("YlOrRd", avg_pm25, n = 5)(avg_pm25),
    weight = 0.5,
    opacity = 1,
    color = "white",
    fillOpacity = 0.7,
    label = ~paste("ZIP:", zcta, "<br>PM2.5:", round(avg_pm25, 2)),
    highlightOptions = highlightOptions(
      weight = 2,
      color = "#666",
      fillOpacity = 0.9,
      bringToFront = TRUE
    )
  ) %>%
  addLegend(
    pal = colorQuantile("YlOrRd", ces_map_data$avg_pm25, n = 5),
    values = ~avg_pm25,
    opacity = 0.7,
    title = "PM2.5",
    position = "bottomright"
  )

###############################################################################
###############################################################################
# -----------------------------------
# Other maps
# -----------------------------------
###############################################################################
###############################################################################

# -----------------------------------
# Prepare joined shape data with ACS + CES metrics
# -----------------------------------

# Ensure zcta columns are characters for merging
ces_full <- ces_full %>%
  mutate(zcta = as.character(zcta))

zcta_shapes_ca <- zcta_shapes_ca %>%
  mutate(zcta = as.character(ZCTA5CE20))

# Join shapefile data with ACS and CES summary data
income_map_data <- left_join(zcta_shapes_ca, ces_full, by = "zcta")

# -----------------------------------
# Static map: Median household income
# -----------------------------------
ggplot(data = income_map_data) +
  geom_sf(aes(fill = income), color = NA) +
  scale_fill_viridis_c(option = "C", na.value = "grey90", name = "Median Income") +
  labs(
    title = "Median Household Income by ZIP Code in California",
    subtitle = "Using ACS 2022 Data",
    caption = "ZIPs with missing data shown in grey"
  ) +
  theme_minimal()

# -----------------------------------
# Static map: Poverty percentage
# -----------------------------------
ggplot(data = income_map_data) +
  geom_sf(aes(fill = pct_poverty), color = NA) +
  scale_fill_viridis_c(option = "C", na.value = "grey90", name = "% Poverty") +
  labs(
    title = "Poverty Percentage by ZIP Code in California",
    subtitle = "Using ACS 2022 Data",
    caption = "ZIPs with missing data shown in grey"
  ) +
  theme_minimal()

# -----------------------------------
# Interactive leaflet map: Income
# -----------------------------------
leaflet(data = income_map_data) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(
    fillColor = ~colorQuantile("YlGnBu", income, n = 5)(income),
    weight = 0.5,
    opacity = 1,
    color = "white",
    fillOpacity = 0.7,
    label = ~paste("ZIP:", zcta, "<br>Income: $", round(income, 0)),
    highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.9, bringToFront = TRUE)
  ) %>%
  addLegend(
    pal = colorQuantile("YlGnBu", income_map_data$income, n = 5),
    values = ~income,
    opacity = 0.7,
    title = "Median Income",
    position = "bottomright"
  )

# -----------------------------------
# Interactive leaflet map: Poverty
# -----------------------------------
leaflet(data = income_map_data) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(
    fillColor = ~colorQuantile("YlOrRd", pct_poverty, n = 5)(pct_poverty),
    weight = 0.5,
    opacity = 1,
    color = "white",
    fillOpacity = 0.7,
    label = ~paste("ZIP:", zcta, "<br>Poverty %:", round(pct_poverty, 1)),
    highlightOptions = highlightOptions(weight = 2, color = "#666", fillOpacity = 0.9, bringToFront = TRUE)
  ) %>%
  addLegend(
    pal = colorQuantile("YlOrRd", income_map_data$pct_poverty, n = 5),
    values = ~pct_poverty,
    opacity = 0.7,
    title = "% Poverty",
    position = "bottomright"
  )

# -----------------------------------
# Define top 50% pollution (median split example)
# -----------------------------------
median_pm25 <- median(income_map_data$avg_pm25, na.rm = TRUE)
income_map_data <- income_map_data %>%
  mutate(top50_pm25 = avg_pm25 > median_pm25)

# Find the median income value (for reference)
median_income <- median(income_map_data$income, na.rm = TRUE)

# -----------------------------------
# GAM model: Poverty ~ smooth(PM2.5)
# -----------------------------------
gam_model <- gam(pct_poverty ~ s(avg_pm25), data = income_map_data)
summary(gam_model)

plot(gam_model, 
     main = "GAM: Relationship between PM2.5 and Poverty",
     xlab = "Average PM2.5",
     ylab = "Poverty (%)",
     shade = TRUE,
     col = "blue")

print(paste("Deviance explained by the model:", round(summary(gam_model)$dev.expl * 100, 1), "%"))

# -----------------------------------
# GAM model: Income ~ smooth(PM2.5)
# -----------------------------------
gam_model_income <- gam(income ~ s(avg_pm25), data = income_map_data)
summary(gam_model_income)

plot(gam_model_income, 
     main = "GAM: Relationship between PM2.5 and Income",
     xlab = "Average PM2.5",
     ylab = "Median Income ($)",
     shade = TRUE,
     col = "blue")

print(paste("Deviance explained by the model:", round(summary(gam_model_income)$dev.expl * 100, 1), "%"))

# -----------------------------------
# Define top 10% pollution
# -----------------------------------
quantile_pm25 <- quantile(income_map_data$avg_pm25, 0.90, na.rm = TRUE)
income_map_data <- income_map_data %>%
  mutate(top10_pm25 = avg_pm25 > quantile_pm25)

# Define bottom 25% income
quantile_income <- quantile(income_map_data$income, 0.25, na.rm = TRUE)
income_map_data <- income_map_data %>%
  mutate(bottom25_income = income < quantile_income)

# Map: Overlap of top 10% pollution & bottom 25% income
ggplot(data = income_map_data) +
  geom_sf(aes(fill = top10_pm25 & bottom25_income), color = NA) +
  scale_fill_manual(values = c("FALSE" = "grey90", "TRUE" = "red"), name = "Overlap") +
  labs(
    title = "ZIP Codes in Top 10% Pollution & Bottom 25% Income",
    subtitle = "California ZIP Code Tabulation Areas",
    caption = "Highlighted areas show high pollution and low income"
  ) +
  theme_minimal()

# Define top 25% poverty
quantile_poverty <- quantile(income_map_data$pct_poverty, 0.75, na.rm = TRUE)
income_map_data <- income_map_data %>%
  mutate(top25_poverty = pct_poverty > quantile_poverty)

# Map: Overlap of top 10% pollution & top 25% poverty
ggplot(data = income_map_data) +
  geom_sf(aes(fill = top10_pm25 & top25_poverty), color = NA) +
  scale_fill_manual(values = c("FALSE" = "grey90", "TRUE" = "red"), name = "Overlap") +
  labs(
    title = "ZIP Codes in Top 10% Pollution & Top 25% Poverty",
    subtitle = "California ZIP Code Tabulation Areas",
    caption = "Highlighted areas show high pollution and high poverty"
  ) +
  theme_minimal()

###############################################################################
###############################################################################
# -----------------------------------
###############################################################################
# --------------------- Valley Fever Part - Clean Rewrite ---------------------
###############################################################################

# ---------------------------------------
# Read and clean Valley Fever data
# ---------------------------------------
valley_fever <- read_csv("ValleyFever_2022_Cases_Rates_Clean.csv")

valley_fever <- valley_fever %>%
  mutate(
    County = str_to_title(County),
    County = str_trim(County),
    Rate_per_100k = as.numeric(Rate_per_100k),
    Cases = as.numeric(Cases)
  )

# ---------------------------------------
# Load and prepare county shapefile
# ---------------------------------------
ca_counties <- counties(cb = TRUE, year = 2020) %>%
  filter(STATEFP == "06") %>%
  mutate(NAME = str_trim(NAME))

ca_counties <- left_join(ca_counties, valley_fever, by = c("NAME" = "County"))

# ---------------------------------------
# Static ggplot map of Valley Fever cases
# ---------------------------------------
ggplot(data = ca_counties) +
  geom_sf(aes(fill = Cases), color = "white") +
  geom_sf_text(aes(label = NAME), size = 2, color = "black", check_overlap = TRUE) +
  scale_fill_viridis_c(option = "C", na.value = "grey90", name = "Cases (2022)") +
  labs(
    title = "Valley Fever Cases by County (2022)",
    subtitle = "California Department of Public Health",
    caption = "Source: CDPH 2022"
  ) +
  theme_minimal()

# ---------------------------------------
# Prepare for interactive leaflet map
# ---------------------------------------
ca_counties <- st_transform(ca_counties, crs = 4326)

pal <- colorNumeric("YlOrRd", domain = ca_counties$Cases, na.color = "transparent")

leaflet(data = ca_counties) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(
    fillColor = ~pal(Cases),
    weight = 1,
    opacity = 1,
    color = "white",
    dashArray = "3",
    fillOpacity = 0.7,
    label = ~paste0(NAME, "<br>Cases: ", Cases, "<br>Rate per 100k: ", round(Rate_per_100k, 1)),
    highlightOptions = highlightOptions(weight = 2, color = "#666", dashArray = "", fillOpacity = 0.9, bringToFront = TRUE)
  ) %>%
  addLegend(
    pal = pal,
    values = ~Cases,
    opacity = 0.7,
    title = "Valley Fever Cases (2022)",
    position = "bottomright"
  )

###############################################################################
# ------------------ ZIP-level high risk + VF overlap ------------------------
###############################################################################

# ---------------------------------------
# ---------------------------------------
# Step 1: Flag high-risk ZIP codes
# ---------------------------------------
pm25_cutoff <- quantile(ces_full$avg_pm25, 0.90, na.rm = TRUE)
income_cutoff <- quantile(ces_full$income, 0.25, na.rm = TRUE)

ces_full <- ces_full %>%
  mutate(
    high_pm25 = avg_pm25 >= pm25_cutoff,
    low_income = income <= income_cutoff,
    high_risk_zip = high_pm25 & low_income
  )

# ---------------------------------------
# Step 2: Flag high VF counties
# ---------------------------------------
vf_rate_cutoff <- quantile(valley_fever$Rate_per_100k, 0.90, na.rm = TRUE)

valley_fever <- valley_fever %>%
  mutate(high_vf = Rate_per_100k >= vf_rate_cutoff)

# ---------------------------------------
# Step 3: Merge county info into ZIP-level data
# ---------------------------------------
zip_county_df <- read_csv("Zip_Code_County_Clean_Final.csv") %>%
  mutate(zip = as.character(zip), county = str_to_title(county))

ces_full <- ces_full %>%
  left_join(zip_county_df, by = c("zcta" = "zip"))

colnames(ces_full)

ces_full <- ces_full %>%
  left_join(valley_fever %>% select(County, high_vf, Rate_per_100k, Cases), 
            by = c("county" = "County")) %>%
  mutate(overlap_high = high_risk_zip & high_vf)

# ---------------------------------------
# Step 4: Join spatial ZIP shape data for mapping
# ---------------------------------------
ces_map_data <- left_join(zcta_shapes_ca, ces_full, by = c("ZCTA5CE20" = "zcta"))

# ---------------------------------------
# Step 5: Create overlap map
# ---------------------------------------
ggplot(data = ces_map_data) +
  geom_sf(aes(fill = overlap_high), color = NA) +
  scale_fill_manual(
    values = c("TRUE" = "red", "FALSE" = "grey80"),
    name = "High PM2.5 & Low Income ZIPs\nin High VF Counties",
    labels = c("No", "Yes")
  ) +
  labs(
    title = "ZIP Codes with High Pollution & Poverty Overlapping High Valley Fever Counties",
    subtitle = "Top 10% PM2.5, Bottom 25% Income, and Top 10% VF Rate",
    caption = "Data: CalEnviroScreen 4.0, ACS 2022, CDPH 2022"
  ) +
  theme_minimal()

# ---------------------------------------
# Step 6: Statistical tests on overlap
# ---------------------------------------
table(ces_full$high_risk_zip, ces_full$overlap_high)
chisq.test(table(ces_full$high_risk_zip, ces_full$overlap_high))
phi(table(ces_full$high_risk_zip, ces_full$overlap_high))

# ---------------------------------------
# Step 7: Proportion plot
# ---------------------------------------
prop_df <- ces_full %>%
  mutate(VF_Label = ifelse(high_vf, "High VF County", "Not High VF County")) %>%
  group_by(VF_Label) %>%
  summarize(prop_high_risk = mean(high_risk_zip, na.rm = TRUE), n = n(), .groups = "drop")

ggplot(prop_df, aes(x = VF_Label, y = prop_high_risk)) +
  geom_point(size = 5, color = "red") +
  geom_segment(aes(xend = VF_Label, yend = 0), color = "grey50") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Proportion of High Risk ZIP Codes by VF County Status", x = "", y = "Proportion of High Risk ZIPs") +
  theme_minimal(base_size = 14)

###############################################################################
# ------------------- County percentile ranks and final join -----------------
###############################################################################

vf_percentiles <- ca_counties %>%
  st_drop_geometry() %>%
  mutate(
    Rate_per_100k = as.numeric(Rate_per_100k),
    percentile_rank = percent_rank(Rate_per_100k) * 100
  ) %>%
  select(NAME, Cases, Rate_per_100k, percentile_rank) %>%
  arrange(desc(Rate_per_100k)) %>%
  mutate(
    percentile_group = case_when(
      percentile_rank >= 90 ~ "Top 10%",
      percentile_rank >= 75 ~ "Top 25%",
      percentile_rank >= 50 ~ "Top 50%",
      TRUE ~ "Bottom 50%"
    )
  )

vf_ranked <- ca_counties %>%
  st_drop_geometry() %>%
  mutate(
    Rate_per_100k = as.numeric(Rate_per_100k),
    percentile_rank = percent_rank(Rate_per_100k)
  ) %>%
  select(NAME, Rate_per_100k, percentile_rank)

valley_fever_clean <- valley_fever %>%
  select(County, Year, Cases, Rate_per_100k) %>%
  left_join(vf_ranked, by = c("County" = "NAME"))

vf_top25_cutoff <- 0.75  # Change this value as needed

valley_fever_clean <- valley_fever_clean %>%
  mutate(top25_vf = percentile_rank >= vf_top25_cutoff)

ces_full <- ces_full %>%
  select(-starts_with("top25_vf"), -starts_with("Rate_per_100k"), -starts_with("percentile_rank"))
colnames(valley_fever_clean)
ces_full <- ces_full %>%
  left_join(valley_fever_clean %>% select(County, top25_vf, Rate_per_100k.y), 
            by = c("county" = "County")) %>%
  mutate(overlap_final = high_risk_zip & top25_vf)

ces_map_data <- left_join(zcta_shapes_ca, ces_full, by = c("ZCTA5CE20" = "zcta"))

ggplot(data = ces_map_data) +
  geom_sf(aes(fill = overlap_final), color = NA) +
  scale_fill_manual(
    values = c("FALSE" = "grey90", "TRUE" = "red"),
    name = "Overlap Status",
    labels = c("No", "Yes")
  ) +
  labs(
    title = "ZIP Codes with High Pollution & Poverty in Top 25% VF Counties",
    subtitle = "Top 10% PM2.5, Bottom 25% Income, and Top 25% Valley Fever Rate Counties",
    caption = "Data: CalEnviroScreen 4.0, ACS 2022, CDPH 2022"
  ) +
  theme_minimal()



# First, compute percentile ranks for counties
vf_percentiles <- ca_counties %>%
  st_drop_geometry() %>%
  mutate(
    Rate_per_100k = as.numeric(Rate_per_100k),
    percentile_rank = percent_rank(Rate_per_100k) * 100,
    percentile_group = case_when(
      percentile_rank >= 90 ~ "Top 10%",
      percentile_rank >= 75 ~ "Top 25%",
      percentile_rank >= 50 ~ "Top 50%",
      TRUE ~ "Bottom 50%"
    )
  ) %>%
  select(NAME, Rate_per_100k, percentile_rank, percentile_group)

# View to confirm
View(vf_percentiles)


valley_fever_clean <- valley_fever %>%
  select(County, Rate_per_100k, Cases) %>%
  left_join(vf_percentiles, by = c("County" = "NAME"))

ces_full <- ces_full %>%
  left_join(valley_fever_clean %>% select(County, percentile_group), 
            by = c("county" = "County"))

ces_map_data <- left_join(zcta_shapes_ca, ces_full, by = c("ZCTA5CE20" = "zcta"))

ggplot(data = ces_map_data) +
  geom_sf(aes(fill = percentile_group), color = NA) +
  scale_fill_manual(
    values = c("Top 10%" = "red", 
               "Top 25%" = "orange", 
               "Top 50%" = "yellow", 
               "Bottom 50%" = "grey80"),
    name = "VF County Group"
  ) +
  labs(
    title = "ZIP Codes by Valley Fever County Percentile Group",
    subtitle = "Based on County VF Rate Percentile",
    caption = "Data: CDPH 2022, CalEnviroScreen 4.0, ACS 2022"
  ) +
  theme_minimal()


vf_percentiles <- ca_counties %>%
  st_drop_geometry() %>%
  mutate(
    Rate_per_100k = as.numeric(Rate_per_100k),
    percentile_rank = percent_rank(Rate_per_100k) * 100,
    percentile_group = case_when(
      percentile_rank >= 90 ~ "Top 10%",
      percentile_rank >= 75 ~ "Top 25%",
      percentile_rank >= 50 ~ "Top 50%",
      TRUE ~ "Bottom 50%"
    )
  ) %>%
  select(NAME, Rate_per_100k, percentile_rank, percentile_group)

valley_fever_clean <- valley_fever %>%
  select(County, Rate_per_100k, Cases) %>%
  left_join(vf_percentiles, by = c("County" = "NAME"))

# Add county column first if needed
zip_county_df <- read_csv("Zip_Code_County_Clean_Final.csv") %>%
  mutate(zip = as.character(zip), county = str_to_title(county))

ces_full <- ces_full %>%
  left_join(zip_county_df, by = c("zcta" = "zip"))
colnames(ces_full)
# Merge percentile group
ces_full <- ces_full %>%
  left_join(valley_fever_clean %>% select(County, percentile_group), by = c("county.y" = "County"))

pm25_cutoff <- quantile(ces_full$avg_pm25, 0.90, na.rm = TRUE)
income_cutoff <- quantile(ces_full$income, 0.25, na.rm = TRUE)

ces_full <- ces_full %>%
  mutate(
    high_pm25 = avg_pm25 >= pm25_cutoff,
    low_income = income <= income_cutoff,
    high_risk_zip = high_pm25 & low_income
  )

ces_map_data <- left_join(zcta_shapes_ca, ces_full, by = c("ZCTA5CE20" = "zcta"))

colnames(ces_map_data)

ces_map_data <- ces_map_data %>%
  mutate(plot_group = ifelse(high_risk_zip, percentile_group.y, "Other"))

ggplot(data = ces_map_data) +
  geom_sf(aes(fill = plot_group), color = NA) +
  scale_fill_manual(
    values = c(
      "Top 10%" = "red",
      "Top 25%" = "orange",
      "Top 50%" = "yellow",
      "Bottom 50%" = "grey60",
      "Other" = "grey90"
    ),
    name = "VF County Group"
  ) +
  coord_sf(
    xlim = c(-124, -114),  # longitude limits (adjust as needed)
    ylim = c(32, 39.5)       # latitude limits (adjust as needed)
  ) +
  labs(
    title = "ZIP Codes by VF County Percentile Group",
    subtitle = "Highlighting ZIPs in Top 10% PM2.5 & Bottom 25% Income",
    caption = "Data: CalEnviroScreen 4.0, ACS 2022, CDPH 2022"
  ) +
  theme_minimal()

###############################################################################
# -------------------------- Valley Fever Data Setup --------------------------
###############################################################################

# Read Valley Fever dataset
valley_fever <- read_csv("ValleyFever_2022_Cases_Rates_Clean.csv") %>%
  mutate(
    County = str_to_title(County),
    County = str_trim(County),
    Rate_per_100k = as.numeric(Rate_per_100k),
    Cases = as.numeric(Cases)
  )

# Load county shapefile and join
ca_counties <- counties(cb = TRUE, year = 2020) %>%
  filter(STATEFP == "06") %>%
  mutate(NAME = str_trim(NAME))

ca_counties <- left_join(ca_counties, valley_fever, by = c("NAME" = "County"))

# View map
ggplot(data = ca_counties) +
  geom_sf(aes(fill = Cases), color = "white") +
  geom_sf_text(aes(label = NAME), size = 2, color = "black", check_overlap = TRUE) +
  scale_fill_viridis_c(option = "C", na.value = "grey90", name = "Cases (2022)") +
  labs(
    title = "Valley Fever Cases by County (2022)",
    subtitle = "California Department of Public Health",
    caption = "Source: CDPH 2022"
  ) +
  theme_minimal()

###############################################################################
# ---------------------- ZIP-level High Risk Flagging -------------------------
###############################################################################

# Create cutoff for PM2.5 and poverty
pm25_cutoff <- quantile(ces_full$avg_pm25, 0.90, na.rm = TRUE)
poverty_cutoff <- quantile(ces_full$pct_poverty, 0.75, na.rm = TRUE)

ces_full <- ces_full %>%
  mutate(
    high_pm25 = avg_pm25 >= pm25_cutoff,
    high_poverty = pct_poverty >= poverty_cutoff,
    high_risk_zip = high_pm25 & high_poverty
  )

###############################################################################
# ---------------------- VF County Percentile Groups --------------------------
###############################################################################

# Create percentile ranks and groups for counties
vf_percentiles <- ca_counties %>%
  st_drop_geometry() %>%
  mutate(
    Rate_per_100k = as.numeric(Rate_per_100k),
    percentile_rank = percent_rank(Rate_per_100k) * 100,
    percentile_group = case_when(
      percentile_rank >= 90 ~ "Top 10%",
      percentile_rank >= 75 ~ "Top 25%",
      percentile_rank >= 50 ~ "Top 50%",
      TRUE ~ "Bottom 50%"
    )
  ) %>%
  select(NAME, Rate_per_100k, percentile_rank, percentile_group)

# Join percentiles to Valley Fever data
valley_fever_clean <- valley_fever %>%
  select(County, Rate_per_100k, Cases) %>%
  left_join(vf_percentiles, by = c("County" = "NAME"))

# Add county column to ZIP-level data if not present
zip_county_df <- read_csv("Zip_Code_County_Clean_Final.csv") %>%
  mutate(zip = as.character(zip), county = str_to_title(county))

ces_full <- ces_full %>%
  left_join(zip_county_df, by = c("zcta" = "zip"))

# Join percentile group to ZIP-level data
ces_full <- ces_full %>%
  left_join(valley_fever_clean %>% select(County, percentile_group), 
            by = c("county.x" = "County"))

###############################################################################
# ---------------------- Map: Highlighting ZIP Codes --------------------------
###############################################################################

# Join to shapes
ces_map_data <- left_join(zcta_shapes_ca, ces_full, by = c("ZCTA5CE20" = "zcta"))

# Create final plot group column
ces_map_data <- ces_map_data %>%
  mutate(plot_group = ifelse(high_risk_zip, percentile_group, "Other"))

ggplot(data = ces_map_data) +
  geom_sf(aes(fill = plot_group), color = NA) +
  scale_fill_manual(
    values = c(
      "Top 10%" = "red",
      "Top 25%" = "orange",
      "Top 50%" = "yellow",
      "Bottom 50%" = "grey60",
      "Other" = "grey90"
    ),
    name = "VF County Group"
  ) +
  coord_sf(
    xlim = c(-124, -114),  # longitude limits (adjust as needed)
    ylim = c(32, 39.5)       # latitude limits (adjust as needed)
  ) +
  labs(
    title = "ZIP Codes by VF County Percentile Group",
    subtitle = "Highlighting ZIPs in Top 10% PM2.5 & Top 25% Poverty",
    caption = "Data: CalEnviroScreen 4.0, ACS 2022, CDPH 2022"
  ) +
  theme_minimal()


###############################################################################
# --------------------- Statistical Tests on Overlap --------------------------
###############################################################################

table(ces_full$high_risk_zip, ces_full$percentile_group)

# You can run more detailed stats as needed (e.g., chi-square, phi coefficient)
