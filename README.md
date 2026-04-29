# ACItools

R package to calculate species-level and community-level Adaptive Capacity Index (ACI) based on species traits and abundance data.

## Installation

```r
install.packages("remotes")
remotes::install_github("LuisEnriqueAngelesGonzalez/ACItools")
```

## Workflow

The package follows four main steps:

1. Prepare and standardise trait data  
2. Select variables using correlation and PCA  
3. Calculate species-level ACI  
4. Calculate community-level ACI  

## Example

```r
library(ACItools)

# Load trait data
df <- read.csv(
  system.file("extdata", "ACI_species.csv", package = "ACItools")
)

# Prepare data
data_scaled <- prepare_trait_data(df)

# Select variables
selected <- select_aci_variables(data_scaled)

# Define ACI variables and direction
aci_variables <- c("K", "M", "RM", "GT", "AOM", "Latitude", "Longitude")
aci_signs <- c(1, 1, 1, -1, -1, 1, 1)

# Calculate species-level ACI
Species_ACI <- calculate_species_aci(
  data_scaled,
  aci_variables,
  aci_signs
)

# Load abundance data
abundance_df <- read.csv(
  system.file("extdata", "ACI_c.csv", package = "ACItools")
)

# Calculate community-level ACI
comm <- calculate_community_aci(
  abundance_df = abundance_df,
  species_aci = Species_ACI
)

df_comm <- comm$community_transect
```

## Notes

- Positive signs indicate variables that increase ACI  
- Negative signs indicate variables that decrease ACI  
- Trait data must include a column named `Species`  
- Abundance data must include `Species`, `Abundance`, `Site_sample`, and `Locality`
