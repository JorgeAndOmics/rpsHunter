# Visual style

Every rpsHunter figure follows one visual system, so a reader who has learnt one page can read any
other. Its executable form is `workflow/scripts/plots/style.R`; plotting code asks that file
instead of deciding for itself.

## Rules

1. **One colour, one meaning.** A recurring concept has the same colour on every page. Plotting
   code never writes a colour value (a test checks).
2. **Colour-blind safe.** The palette is Okabe and Ito's. Every pair of its hues stays
   distinguishable under deuteranopia, protanopia and tritanopia (smallest CIELAB distance 16.2;
   a test holds it above 15).
3. **Ordered things get ordered colours.** An ordered set takes a dark-to-light ramp, so "stronger"
   reads as darker before the legend is consulted.
4. **Species on rows, always in the same order,** and every configured species keeps its row, empty
   if it has no result.
5. **Pages scale.** A page's height follows the number of species, so 100 species are as readable
   as 5.

## Colour

| Name | Hex |
|---|---|
| blue | `#0072B2` |
| orange | `#E69F00` |
| green | `#009E73` |
| purple | `#CC79A7` |
| vermillion | `#D55E00` |
| sky | `#56B4E9` |
| yellow | `#F0E442` |

Fixed meanings:

- **Evidence tier of a CDD call** (standard, borderline, rescue): a ramp from the palette's blue to
  a pale tint, strongest darkest.
- **How CDD classifies a call** (specific, non-specific, superfamily): the same ramp.
- **Method agreement**: both methods green, CDD only orange, HMMER only sky, not compared grey.
- **Counts on tiles**: a pale-to-blue scale, logarithmic, so that a family with ten loci is visible
  beside one with a thousand. Numbers are printed on the tiles.
- **Domain families** have no fixed meaning (they are yours), so they take the palette's hues in
  config order, always with a legend. Past seven families, later ones take lighter tints.

Neutrals: grey for "not compared" and for the backbone of a locus; near-black ink for text.

## Type and words

- **IBM Plex Sans** everywhere, installed with the conda environment and embedded in the PDFs, so
  every machine renders the same figure.
- Titles are bold, left-aligned, sentence case. The subtitle says what the page shows and how to
  read it.
- Species binomials are italic.
- Identifiers are shown as words (`cdd_only` as "CDD only"); underscores become spaces.
- Counts are whole numbers with thousands separators.

## Species order

The order of `species:` in the config, or, when `figures.species_tree` names a Newick file, the
order of the tree's tips (matched by species key or display name; species the tree lacks go last).
The tree gives the order only; it is not drawn.

## Output

Four vector PDFs under `results/plots/`, described in [outputs.md](outputs.md#figures). The first
three open with a key to the colours. A PDF has one sheet size, set by its largest page; the pages
of one file are sized alike.

## Adding a page

1. Prepare its counts in `plots/tables.R` as a function of the run's tables, with a test.
2. Draw it in `plots/pages.R` as `page_<name>(run, layout)`, returning `page(plot, size)`. Use the
   scales and the theme of `style.R`; write no colour.
3. Add it to `pdf_pages()` in `figures.R` and to the synthetic-run test, which draws every page at
   5, 30 and 100 species and for a run with no domain.
