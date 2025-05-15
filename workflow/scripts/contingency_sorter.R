# =============================================================================
# DEPENDENCIES
# =============================================================================
options(warn = -1)
suppressMessages({
  library(arrow)
  library(tidyverse)
  library(ggsci)
  library(plotly)
})

# =============================================================================
# COMMAND LINE ARGUMENTS
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

args.data <- args[1]
args.species <- args[2]
args.output_plot_folder <- args[3]
args.output_table_folder <- args[4]

# =============================================================================
# MESSAGE
# =============================================================================
print("Parsing domain data...")

# =============================================================================
# DATA IMPORT
# =============================================================================
data <- arrow::read_parquet(args.data)

# Species
full.species <- readLines(args.species)
full.species <- gsub("_", " ", full.species)

# Clean empty strings and data types
data.clean <- data %>%
  dplyr::rename(Species = File_Name, Domain = Short_name, Bitscore = Bit_score) %>%
  filter(nzchar(Domain)) %>%
  mutate(
    Species = gsub("_", " ", Species),
    Bitscore = as.numeric(Bitscore)
  )

# =============================================================================
# DATA PROCESSING
# =============================================================================

# Both N and C hold the same degree of incompleteness plot-wise, so we replace
# them with something more descriptive
data.group <- data.clean %>%
  group_by(Species, Domain) %>%
  summarise(
    Bitscore = mean(Bitscore),
    Incomplete = case_when(
      any(Incomplete == "-") ~ "W",
      any(Incomplete == "NC") ~ "NC",
      any(Incomplete == "N") & any(Incomplete == "C") ~ "B",
      all(Incomplete == "N") ~ "N",
      all(Incomplete == "C") ~ "C"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Incomplete.desc = case_when(
      Incomplete == "W"  ~ "Complete",
      Incomplete == "N"  ~ "Truncated",
      Incomplete == "C"  ~ "Truncated",
      Incomplete == "B"  ~ "Truncated",
      Incomplete == "NC" ~ "Bitruncated"
    )
  )

# Extract unique domains from data
unique.domains <- data.group %>%
  distinct(Domain) %>%
  pull(Domain)
unique.domains <- unique.domains[nzchar(unique.domains)]

# Extract non-empty species from user-provided file
full.species <- full.species[nzchar(full.species)]

# =============================================================================
# MISSING DATA IMPUTATION (CARTESIAN PRODUCT)
# =============================================================================

# Existing data
existing.data <- data.group

# Species with no hits
missing.species <- full.species[!full.species %in% intersect(full.species, existing.data$Species)]

# Generate cartesian product of all species vs. domains
missing.data <- expand.grid(
  Species = full.species,
  Domain = unique.domains,
  Incomplete = "Z",
  Incomplete.desc = "No Hit",
  Bitscore = 0
)

# Combine existing + missing
plot.data <- bind_rows(existing.data, missing.data)

# =============================================================================
# RANKING AND SLICE
# =============================================================================

domain.ranking <- c("Complete" = 1,
                    "Truncated" = 2,
                    "Bitruncated" = 3,
                    "No Hit" = 4)

plot.data.ranked <- plot.data %>%
  mutate(
    rank = domain.ranking[Incomplete.desc],
    Domain = factor(Domain, levels = rev(sort(unique(Domain))))
  )

# =============================================================================
# TILE PLOT (ggplot2)
# =============================================================================

completion.levels <- c("Complete", "Truncated", "Bitruncated", "No Hit")
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
  geom_text(
    aes(label = ifelse(Incomplete %in% c("N", "C", "B"), Incomplete, "")),
    alpha = 1, color = "black", size = 2
  ) +
  scale_fill_manual(values = custom_colors) +
  labs(fill = "Domain Completion") +
  theme_minimal() +
  guides(
    alpha = FALSE,
    fill = guide_legend(
      title.theme = element_text(face = "bold"),
      label.theme = element_text(face = "bold")
    )
  ) +
  theme(
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 7),
    axis.text.x = element_text(
      face = "bold.italic",
      size = 7,
      angle = 70,
      hjust = 1
    )
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
