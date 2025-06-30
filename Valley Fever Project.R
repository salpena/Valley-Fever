# Clear all objects from the workspace
rm(list = ls())

# Install required packages
install.packages("tidycensus")
install.packages("zctaCrosswalk")
install.packages("readxl")

# Load libraries
library(tidyverse)    # Collection of R packages for data science
library(readr)        # For reading CSV files
library(zctaCrosswalk)# For working with ZIP Code Tabulation Areas (ZCTAs)
library(tidycensus)   # To access Census data via API
library(dplyr)        # For data manipulation
library(readxl)       # For reading Excel files

# Check tidycensus package (just shows its content)
head(tidycensus)

# Set Census API key (replace with your own key if needed)
census_api_key("94dc655738967530ffcdc4dc68b34b420b902d6a", install = TRUE)

# Reload environment so the key is active immediately
readRenviron("~/.Renviron")  

# Re-load these libraries to ensure they're active
library(tidycensus)
library(dplyr)
library(readr)

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
View(CalEnviroScreen_4_0_Results)  # View entire dataset

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
head(ces_clean$tract)

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

# Load ggplot2 for plotting
library(ggplot2)

# Create bar plot of top 20 ZIP codes by PM2.5
ces_full %>%
  arrange(desc(avg_pm25)) %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = reorder(zcta, avg_pm25), y = avg_pm25)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(title = "Top 20 Most Polluted ZIP Codes (PM2.5)",
       x = "ZIP Code", y = "Average PM2.5")

# Install and load spatial packages
install.packages("tigris")
install.packages("sf")
library(tigris)
library(sf)

# Enable caching for shapefiles
options(tigris_use_cache = TRUE)

# Download ZCTA shapefiles
zcta_shapes <- zctas(cb = TRUE, year = 2020)
colnames(zcta_shapes)

# Optionally filter for ZIP codes starting with "9" (California)
zcta_shapes <- zctas(cb = TRUE, year = 2020) %>%
  filter(startsWith(ZCTA5CE200, "9"))

# Filter shapes to ZIP codes present in CES data
zcta_shapes_ca <- zcta_shapes %>%
  filter(ZCTA5CE20 %in% ces_zip_summary$zcta)

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

# Install and load leaflet for interactive maps
install.packages("leaflet")
library(leaflet)

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
