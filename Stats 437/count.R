# lyme_lab_simple.R
# Short lab: Poisson vs Negative Binomial for Lyme disease counts

# --- Load libraries --------------------------------------------------------
library(tidyverse) # data wrangling & plotting
library(sf)          # spatial geometry handling (used for area / centroids)
library(usmap)       # county/state shapes
library(MASS)        # glm.nb (Negative Binomial)
library(AER)         # dispersiontest
library(broom)       # tidy model output

# --- Make temperature data--------------------------------------------------

state_temp <- tribble(
  ~stname,          ~mean_temp_F,
  "Alabama",        63.5,
  "Alaska",         28.1,
  "Arizona",        61.1,
  "Arkansas",       61.1,
  "California",     59.9,
  "Colorado",       46.4,
  "Connecticut",    49.9,
  "Delaware",       56.1,
  "Florida",        71.5,
  "Georgia",        64.3,
  "Hawaii",         70.2,
  "Idaho",          45.1,
  "Illinois",       52.7,
  "Indiana",        52.5,
  "Iowa",           48.7,
  "Kansas",         55.2,
  "Kentucky",       56.5,
  "Louisiana",      67.2,
  "Maine",          41.7,
  "Maryland",       55.4,
  "Massachusetts",  48.8,
  "Michigan",       45.3,
  "Minnesota",      42.2,
  "Mississippi",    64.2,
  "Missouri",       55.5,
  "Montana",        43.3,
  "Nebraska",       49.6,
  "Nevada",         51.1,
  "New Hampshire",  44.6,
  "New Jersey",     53.6,
  "New Mexico",     53.8,
  "New York",       46.1,
  "North Carolina", 59.9,
  "North Dakota",   41.3,
  "Ohio",           51.7,
  "Oklahoma",       60.6,
  "Oregon",         49.3,
  "Pennsylvania",   49.4,
  "Rhode Island",   50.9,
  "South Carolina", 63.3,
  "South Dakota",   46.1,
  "Tennessee",      58.5,
  "Texas",          65.3,
  "Utah",           49.5,
  "Vermont",        43.2,
  "Virginia",       56.3,
  "Washington",     49.1,
  "West Virginia",  52.7,
  "Wisconsin",      44.3,
  "Wyoming",        42.8
)


# --- Read data and basic transforms ---------------------------------------
# Expect lyme_data.csv to contain Cases2022, TOT_POP, stcode, ctycode, and other county fields.
lyme <- read_csv("lyme_data.csv")

# Create standard FIPS code and a simple case rate
lyme <- lyme %>%
  mutate(
    FIPS = paste0(sprintf("%02d", as.numeric(stcode)),
                  sprintf("%03d", as.numeric(ctycode))),
    case_rate_2022 = Cases2022 / TOT_POP
  )

# --- State-level mean temperature (paste your state_temp tibble here) -----
# Example: state_temp <- tribble( ~stname, ~mean_temp_F, "Maine", 41, ... )
# If you already created state_temp in the slides, this join will use it.
# (Insert your full state_temp tibble above if not already present.)

# --- Get county geometry and merge ---------------------------------------
county_geo <- us_map("counties") %>% st_as_sf()

# Join lyme data to county geometry; keep required variables for modeling.
lyme_geo <- county_geo %>%
  left_join(lyme, by = c("fips" = "FIPS")) %>%
  left_join(state_temp, by = "stname") %>% 
  mutate(
    centroid = st_centroid(geom),
    lat = st_coordinates(centroid)[,2],
    area_km2 = as.numeric(st_area(geom)) / 1e6,
    log_pop_density = log(TOT_POP / area_km2),
    case_rate_2022 = Cases2022 / TOT_POP,
    Northeast = stname %in% c("Massachusetts","Connecticut","New York",
                              "New Jersey","Pennsylvania","Rhode Island",
                              "Vermont","New Hampshire","Maine"
    )
  )  %>% 
  filter(
    !is.na(Cases2022),
    !is.na(TOT_POP),
    !is.na(mean_temp_F),
    !is.na(lat)
  )
# --- Quick EDA plots (students will discuss these in groups) --------------
# 1) Map of incidence rate
usa_data <- us_map("counties") %>% 
  left_join(lyme, by = c("fips" = "FIPS"))

usa_states <- us_map("states")

ggplot(usa_data) + 
  theme_bw() +
  geom_sf(aes(fill = case_rate_2022)) +
  geom_sf(data = usa_states, fill = NA, col = "black", linewidth = 0.35) +
  scale_fill_gradientn(
    colors = c("white", "orange", "magenta", "darkred"),
    labels = scales::label_percent()
  ) +
  labs(
    fill = "Annual\nIncidence\nRate",
    title = "Lyme Disease per-year Incidence Rate (2022)"
  )

# 2) Rate vs mean temperature
ggplot(lyme_geo, aes(x = mean_temp_F, y = case_rate_2022, color = Northeast)) +
  geom_point(alpha = 0.6) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Mean state temperature (F)", y = "Incidence rate", title = "Incidence rate vs mean temperature")

# 3) Count distribution
ggplot(lyme, aes(x = Cases2022)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  labs(title = "Distribution of county Lyme case counts (2022)", x = "Cases (2022)", y = "Count of counties")

# --- Modeling: Poisson with offset (rates) --------------------------------
# We model counts with log(TOT_POP) as an offset to estimate incidence rates.
fit_pois <- glm(
  Cases2022 ~ ststatus + mean_temp_F + log_pop_density + offset(log(TOT_POP)),
  family = poisson(link = "log"),
  data = lyme_geo
)

# Summarize results (students interpret IRRs)
tidy(fit_pois, exponentiate = TRUE, conf.int = TRUE)

# --- Check for overdispersion ---------------------------------------------
# Simple diagnostic: residual deviance / residual df
disp_stat <- deviance(fit_pois) / df.residual(fit_pois)
disp_stat

# Another check: compare mean and variance of counts
mean(lyme$Cases2022, na.rm = TRUE)
var(lyme$Cases2022, na.rm = TRUE)

# Formal test (AER::dispersiontest)
dispersiontest(fit_pois)

# --- Fit Negative Binomial (if overdispersed) -----------------------------
fit_nb <- glm.nb(
  Cases2022 ~ mean_temp_F + log_pop_density + offset(log(TOT_POP)) + stname,
  data = lyme_geo
)

fit_nb_poly <- glm.nb(
  Cases2022 ~ poly(mean_temp_F,2) + log_pop_density + offset(log(TOT_POP)) + stname,
  data = lyme_geo
)


summary(fit_nb)

tidy(fit_nb, exponentiate = TRUE, conf.int = TRUE)

# --- Model comparison: AIC and LRT ----------------------------------------
AIC(fit_pois, fit_nb)

# Likelihood ratio test (NB vs Poisson)
# Note: using anova() between glm and glm.nb is allowed; specify test = "LRT".
# but you'll need to get the p-value by hand.

anova_out = anova(fit_pois, fit_nb, test = "LRT")
1 - pchisq(anova_out$Deviance[2],df = 1)

# --- Predictions & visualization ------------------------------------------
# Predicted rates from the preferred model (use NB if chosen)
lyme_geo <- lyme_geo %>%
  mutate(pred_nb = predict(fit_nb, type = "response"),
         pred_pois = predict(fit_pois, type = "response"))

# Plot observed vs predicted (NB)
ggplot(lyme_geo, aes(x = TOT_POP, y = pred_nb)) +
  geom_point(alpha = 0.4) +
  labs(title = "Predicted counts (Negative Binomial) vs population", x = "Population", y = "Predicted counts")

# Map predicted incidence rates (NB)
lyme_geo <- lyme_geo %>% mutate(pred_rate_nb = pred_nb / TOT_POP)

ggplot(lyme_geo) +
  geom_sf(data = us_map("states") %>% st_as_sf(), fill = "white", color = "white", size = 0.25) +
  geom_sf(aes(fill = pred_rate_nb), color = NA) +
  scale_fill_viridis_c(name = "Predicted rate", option = "B", labels = scales::label_percent(accuracy = 0.01)) +
  labs(title = "Predicted Lyme incidence rate (Negative Binomial)")