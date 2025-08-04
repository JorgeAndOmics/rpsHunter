# =============================================================================
# DEPENDENCIES
# =============================================================================
options(warn = -1)
suppressMessages({
  library(arrow)
  library(tidyverse)
  library(ggsci)
  library(plotly)
  library(GenomicRanges)
  library(plyranges)
})

# =============================================================================
# COMMAND LINE ARGUMENTS
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

args.data <- 'C:/Users/Lympha/Documents/Repositories/rpsHunter/results/tables/domains.parquet'
args.species <- bat_species <- c(
  "Desmodus_rotundus",
  "Miniopterus_schreibersii",
  "Tadarida_brasiliensis",
  "Antrozous_pallidus",
  "Molossus_molossus",
  "Artibeus_lituratus",
  "Eptesicus_fuscus",
  "Myotis_myotis",
  "Eptesicus_nilssonii",
  "Pipistrellus_kuhlii",
  "Rhinolophus_ferrumequinum",
  "Saccopteryx_bilineata",
  "Vespertilio_murinus",
  "Plecotus_auritus",
  "Rhinolophus_hipposideros",
  "Phyllostomus_discolor",
  "Myotis_daubentonii",
  "Myotis_mystacinus",
  "Corynorhinus_townsendii",
  "Hipposideros_larvatus",
  "Rhynchonycteris_naso",
  "Saccopteryx_leptura",
  "Molossus_alvarezi",
  "Glossophaga_mutica",
  "Molossus_nigricans"
)
args.output_plot_folder <- 'C:/Users/Lympha/Desktop/script_test'
args.output_table_folder <- 'C:/Users/Lympha/Desktop/script_test'

# =============================================================================
# MESSAGE
# =============================================================================
print("Parsing domain data...")

# =============================================================================
# DATA IMPORT
# =============================================================================
data <- arrow::read_parquet(args.data)

# Species
full.species <- args.species

# Clean empty strings and data types
data.clean <- data %>%
  filter(nzchar(Domain)) %>%
  mutate(
    Bitscore = as.numeric(Bitscore)
  )

# =============================================================================
# RANGES PROCESSING
# =============================================================================
data.ranges <- GRanges(
  seqnames = data.clean$Chromosome,
  ranges = IRanges(start = data.clean$Start, end = data.clean$End)
  species = data.clean$Species,
  domain = data.clean$Domain,
  bitscore = data.clean$Bitscore,
  evalue = data.clean$Evalue,
  incomplete = data.clean$Incomplete
  pssm = data.clean$PSSM_ID
  superfamily _= data.clean$Superfamily_PSSM_ID
)

# =============================================================================
# FULL CONTINGENCY TABLE & 3D SCATTER PLOT
# =============================================================================

data.grouped.contingency <- data.clean %>%
  group_by(Species, Domain, Incomplete) %>%
  arrange(.by_group = TRUE) %>%
  summarise(Count = n(), .groups = "drop")

incomplete.vec <- c("N", "C", "NC", "-")
domain.vec <- unique(data.clean$Domain)
species.vec <- unique(data.clean$Species)

data.full.cartesian <- expand.grid(
  Species = species.vec,
  Domain = domain.vec,
  Incomplete = incomplete.vec,
  Count = 0
)

# Left-join to get full contingency
data.full.contingency <- data.full.cartesian %>%
  left_join(data.grouped.contingency, by = c("Species", "Domain", "Incomplete")) %>%
  mutate(Count = ifelse(is.na(Count.y), 0, Count.y)) %>%
  select(-c(Count.x, Count.y)) %>%
  mutate(
    Incomplete = case_when(
      Incomplete == "-"  ~ "Complete",
      Incomplete == "NC" ~ "Bitruncated",
      Incomplete == "N"  ~ "N-Truncated",
      Incomplete == "C"  ~ "C-Truncated"
    )
  ) %>%
  arrange(Species, Domain, Incomplete)

# Prepare the 3D data
df.3d <- data.full.contingency %>%
  filter(Count > 0) %>%
  mutate(
    Species = factor(Species),
    Domain = factor(Domain),
    Incomplete = factor(
      Incomplete,
      levels = c("Complete", "N-Truncated", "C-Truncated", "Bitruncated")
    )
  )

# Numeric positions for axes
species_levels <- levels(df.3d$Species)
domain_levels <- levels(df.3d$Domain)
incomplete_levels <- levels(df.3d$Incomplete)

df.3d.hover <- df.3d %>%
  mutate(
    x = as.integer(Species),
    y = as.integer(Domain),
    z = as.integer(Incomplete),
    hover = paste(
      "Species:", Species, "<br>",
      "Domain:", Domain, "<br>",
      "Incomplete:", Incomplete, "<br>",
      "Count:", Count
    )
  )

# Create the 3D scatter plot
scatter.3D.plot <- plot_ly(
  df.3d.hover,
  x = ~x, y = ~y, z = ~z,
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = ~sqrt(Count) * 3,
    color = ~Count,
    colorscale = 'Plasma',
    showscale = TRUE,
    opacity = 0.8
  ),
  text = ~hover,
  hoverinfo = "text"
) %>%
  layout(
    scene = list(
      xaxis = list(
        title = "",
        tickvals = seq_along(species_levels),
        ticktext = species_levels
      ),
      yaxis = list(
        title = "Domain",
        tickvals = seq_along(domain_levels),
        ticktext = domain_levels
      ),
      zaxis = list(
        title = "Completeness",
        tickvals = seq_along(incomplete_levels),
        ticktext = incomplete_levels
      )
    ),
    title = "3D Scatter Plot of Species × Domain × Incomplete × Count"
  )

# =============================================================================
# EXPORTS (Save plots and tables)
# =============================================================================

ggsave(
  file.path(args.output_plot_folder, "tile_plot.png"),
  tile.plot,
  width = 15,
  height = 10,
  dpi = 300
)

htmlwidgets::saveWidget(
  scatter.3D.plot,
  file.path(args.output_plot_folder, "scatter_3D_plot.html")
)

write_csv2(data.full.contingency, file.path(args.output_table_folder, "contingency_table.csv"))
