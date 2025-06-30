rm(list = ls())
install.packages("tidycensus")
install.packages("zctaCrosswalk")
install.packages("readxl")
library(tidyverse)
library(zctaCrosswalk)
library(tidycensus)
library(dplyr)
library(readxl)
head(tidycensus)
census_api_key("94dc655738967530ffcdc4dc68b34b420b902d6a", install = TRUE)
readRenviron("~/.Renviron")  # Load it immediately
library(tidycensus)
library(dplyr)
library(readr)

vars <- c(
  income = "B19013_001",       # Median household income
  pop_total = "B01003_001",    # Total population
  white = "B02001_002",        # White alone
  black = "B02001_003",        # Black or African American alone
  hispanic = "B03003_003",     # Hispanic or Latino
  poverty = "B17001_002"       # Below poverty line
)

zcta_data <- get_acs(
  geography = "zcta",
  variables = vars,
  year = 2022,
  survey = "acs5",
  output = "wide"
)

ca_zctas <- get_zctas_by_state("CA")

ca_data <- zcta_data %>%
  filter(GEOID %in% ca_zctas)

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

ca_data_clean <- ca_data_clean %>%
  mutate(
    pct_white = white / total_pop * 100,
    pct_black = black / total_pop * 100,
    pct_hispanic = hispanic / total_pop * 100,
    pct_poverty = poverty / total_pop * 100
  )

# View the first 10 rows in the console
head(ca_data_clean, 10)

# Or use View() in RStudio
View(ca_data_clean)

# Replace the path with your actual file name
ces_data <- read_csv("Valley Fever Project/CalEnviroScreen_4.0_Results.csv")
View(CalEnviroScreen_4_0_Results)

# View the column names
colnames(ces_data)

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

head(ces_clean)

crosswalk <- read_excel("Valley Fever Project/ZIP_TRACT_122022.xlsx")

colnames(crosswalk)
head(crosswalk)

crosswalk_ca <- crosswalk %>%
  filter(USPS_ZIP_PREF_STATE == "CA") %>%
  mutate(
    tract = as.character(TRACT),
    zcta = as.character(ZIP)
  ) %>%
  select(tract, zcta)

# What do the tract values look like?
head(crosswalk_ca$tract)
head(ces_clean$tract)

crosswalk_ca <- crosswalk_ca %>%
  mutate(tract = sprintf("%011s", tract))

ces_clean <- ces_clean %>%
  mutate(tract = sprintf("%011s", tract))

ces_zip <- left_join(crosswalk_ca, ces_clean, by = "tract")

summary(ces_zip$pm25)

head(ces_zip)

ces_zip_summary <- ces_zip %>%
  group_by(zcta) %>%
  summarize(
    avg_pm25 = mean(pm25, na.rm = TRUE),
    avg_ces_score = mean(ces_score, na.rm = TRUE),
    avg_pct_poverty = mean(pct_poverty, na.rm = TRUE),
    avg_pct_hispanic = mean(pct_hispanic, na.rm = TRUE),
    .groups = "drop"
  )

ces_full <- inner_join(ca_data_clean, ces_zip_summary, by = "zcta")

library(ggplot2)

ces_full %>%
  arrange(desc(avg_pm25)) %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = reorder(zcta, avg_pm25), y = avg_pm25)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(title = "Top 20 Most Polluted ZIP Codes (PM2.5)",
       x = "ZIP Code", y = "Average PM2.5")
install.packages("tigris")
install.packages("sf")
library(tigris)
library(sf)



# Download all ZCTAs but immediately filter only to ZIPs in your dataset
options(tigris_use_cache = TRUE)
zcta_shapes <- zctas(cb = TRUE, year = 2020)
colnames(zcta_shapes)

# Download all ZCTAs, but only keep ZIPs that start with "9"
zcta_shapes <- zctas(cb = TRUE, year = 2020) %>%
  filter(startsWith(ZCTA5CE200, "9"))
# Filter just the ZIPs in your dataset
zcta_shapes_ca <- zcta_shapes %>%
  filter(ZCTA5CE20 %in% ces_zip_summary$zcta)

ces_zip_summary <- ces_zip_summary %>%
  mutate(zcta = as.character(zcta))

zcta_shapes_ca <- zcta_shapes_ca %>%
  mutate(zcta = as.character(ZCTA5CE20))

ces_map_data <- left_join(zcta_shapes_ca, ces_zip_summary, by = "zcta")


library(ggplot2)

ggplot(data = ces_map_data) +
  geom_sf(aes(fill = avg_pm25), color = NA) +
  scale_fill_viridis_c(option = "C", na.value = "grey90", name = "PM2.5") +
  labs(
    title = "Air Pollution (PM2.5) by ZIP Code in California",
    subtitle = "Using CalEnviroScreen 4.0 + ZIP Crosswalk",
    caption = "ZIPs with missing data shown in grey"
  ) +
  theme_minimal()

install.packages("leaflet")
library(leaflet)


leaflet(data = ces_map_data) %>%
  addProviderTiles("CartoDB.Positron") %>%  # clean base map
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
