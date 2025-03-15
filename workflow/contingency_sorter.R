
# DEPENDENCIES
options(warn = -1)
suppressMessages({
  library(arrow)
  library(tidyverse)
  library(ggsci)
})


# COMMAND LINE ARGUMENTS
args <- commandArgs(trailingOnly = TRUE)

args.data <- args[1]
# args.data <- file.path("C:/Users/Lympha/Desktop/results-PRDM9-isoformA/tables/domains.parquet")
args.species <- args[2]
# args.species <- file.path("C:/Users/Lympha/Documents/Repositories/rpsHunter/data/config/species.txt")
args.output_folder <- args[3]
# args.output_folder <- file.path("C:/Users/Lympha/Desktop/results-PRDM9-isoformA/plots")


# MESSAGE
print("Generating tile plot for domains")


# DATA IMPORT
# Domains
data <- arrow::read_parquet(args.data)


# Species
full.species <- readLines(args.species)
full.species <- gsub("_", " ", full.species)

## Clean empty strings and data types
data.clean <- data %>%
  dplyr::rename(Species = File_Name, Domain = Short_name, Bitscore = Bit_score) %>%
  filter(nzchar(Domain)) %>%
  mutate(
    Species = gsub("_", " ", Species),
    Bitscore = as.numeric(Bitscore))


# DATA PROCESSING
## Both N and C hold the same degree of incompleteness plot-wise, so we are changing
## the column values for something more descriptive

data.group <- data.clean %>%
  group_by(Species, Domain) %>%
  summarise(
    Bitscore = mean(Bitscore),
    Incomplete = case_when(
      any(Incomplete == "-") ~ "W",
      any(Incomplete == "NC") ~ "NC",
      any(Incomplete == "N") & any(Incomplete == "C") ~ "B",
      all(Incomplete == "N") ~ "N",
      all(Incomplete == "C") ~ "C",
  )
  ) %>%
  mutate(Incomplete.desc = case_when(
    Incomplete == "W" ~ "Complete",
    Incomplete == "N" ~ "Truncated",
    Incomplete == "C" ~ "Truncated",
    Incomplete == "B" ~ "Truncated",
    Incomplete == "NC" ~ "Bitruncated",
  )) %>%
  ungroup()


# Extract unique values for domain (Domain) from data
unique.domains <- data.group %>%
  distinct(Domain) %>%
  pull(Domain)
unique.domains <- unique.domains[nzchar(unique.domains)] # Select non-empty strings

# Extract unique values for species 
full.species <- full.species[nzchar(full.species)] # Select non-empty strings


# MISSING DATA IMPUTATION
## Determine exiting data
existing.data <- data.group

## Determine species with no hits (absent in data)
missing.species <- full.species[!full.species %in% intersect(full.species, existing.data$Species)]

## Generate Cartesian product of missing species and unique domains
missing.data <- expand.grid(
  Species = full.species,  # We use full species to have a grid in the plot
  Domain = unique.domains,
  Incomplete = "Z",
  Incomplete.desc = "No Hit",
  Bitscore = 0
)

## Add missing data to plot data
plot.data <- bind_rows(existing.data, missing.data)


# GENERATE RANKING AND SLICE
domain.ranking <- c("Complete" = 1, 
                    "Truncated" = 2,
                    "Bitruncated" = 3,
                    "No Hit" = 4)

plot.data.ranked <- plot.data %>%
  mutate(rank = domain.ranking[Incomplete.desc],
         Domain = factor(Domain, levels = sort(unique(Domain)))) 


# GENERATE TILE PLOT
## Define factor for legend order
completion.levels <- c("Complete", "Truncated", "Bitruncated", "No Hit")

# Define custom colors using palette and override "Blank" with white
custom_colors <- pal_bmj()(length(completion.levels))
names(custom_colors) <- completion.levels
custom_colors["No Hit"] <- "white"

tile.plot <- ggplot(plot.data.ranked) +
  aes(
    x = Species,
    y = Domain,
    fill = factor(Incomplete.desc, levels = completion.levels),
    alpha = Bitscore
  ) +
  geom_tile(color = "black") +
  geom_text(aes(label = ifelse(Incomplete == "N" | Incomplete == "C" | Incomplete == "B", Incomplete, "")), 
            alpha = 1, color = "black", size = 2) +
  scale_fill_manual(values = custom_colors) +
  labs(fill = "Domain Completion") +
  theme_minimal() +
  guides(alpha = FALSE,
         fill = guide_legend(
           title.theme = element_text(face = "bold"),
           label.theme = element_text(face = "bold"))) +
  theme(
    axis.title.x = element_text(face = "bold", size = 12L),
    axis.title.y = element_text(face = "bold", size = 12L),
    axis.text.y = element_text(face = "bold", size = 7L),
    axis.text.x = element_text(
      face = "bold.italic",
      size = 7L,
      angle = 70L,
      hjust = 1L
    ),
    
  )

ggsave(file.path(args.output_folder, "tile_plot.png"), tile.plot, width = 15, height = 10, dpi = 300)
