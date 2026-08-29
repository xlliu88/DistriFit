# PRN Generator & Distribution Fitting — Shiny App

An interactive tool to fit data into probability distributions.

Supported distributions:

- normal

## Files

```
ui.R        - page layout / inputs / outputs
server.R    - reactive logic: builds the dataset(s), runs distrEstimate(),
              picks the top fits, renders the plots / table / text
global.R    - sources plus Shiny-only
              helpers: default parameters for the "generate data" panel
              (read from R/distribution.names.txt), a density/pmf lookup,
              dynamic histogram bin selection, and the plotting functions
R/          - your original source files, plus the fixes noted below
```

## Issues: 
1. add metaInfo variable
2. change name convention. pass metaInfo for plot information display
3. reformat fitting result
5. sometimes it doesn't produce density plot after fitting

## Running it

Requires R with the packages the original code already depends on:
`shiny`, `dplyr`, `stringr`, `rlang`, `randtests`, `EnvStats`.

```r
install.packages(c("shiny", "dplyr", "stringr", "rlang", "randtests", "EnvStats"))
shiny::runApp("path/to/this/folder")
```

## Testing note

Still no CRAN access in this sandbox, so — same as last round — the
actual Shiny UI/reactivity couldn't be click-tested in a browser. What
*was* verified directly with a real R 4.3 interpreter and small
stand-ins for `dplyr`/`stringr`/`rlang`:

- Simulated the exact namespaced-ID generation flow for `norm` and
  `pois` (and confirmed `binom`'s IDs and `norm`'s IDs are now
  structurally distinct strings, so no collision is possible).
- Simulated a 3-numeric-column upload, fit one column, and rendered all
  three preview panels — two plain histograms, one with the fitted
  curves overlaid on the *same* histogram (single plot, titled with the
  column name) — confirmed visually.
- Re-ran the full 14-distribution generate → fit → plot sweep with the
  restored descending/`p > cutoff` logic and the new `title` parameter;
  every distribution's own generated data now correctly ranks itself
  at or near the top of its own fit (a good sanity check that the
  ranking direction is right again).
- All input/output IDs cross-checked between `ui.R` and `server.R`
  (including the new dynamic `fit_column`, `seriesPlot_*`,
  `seriesStats_*` outputs, which are necessarily unchecked by static
  grep since they're built at runtime — reviewed by hand instead).
- Every file parses without syntax errors under R 4.3.

Not tested: live Shiny reactivity (the dynamic per-column
`output[[...]] <- renderPlot(...)` registration pattern in particular is
correct Shiny usage but depends on the reactive graph behaving as
expected in a real session), `downloadHandler`, and the clipboard JS.
