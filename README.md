# PRN Generator & Distribution Fitting — Shiny App

An interactive Shiny front end for the existing PRN-generation / chi-squared
goodness-of-fit code (`R/utilities.r`, `R/prn_generator.r`, `R/prob_fit.r`,
`R/gofTest.chisq.r`, `R/main.r`).

## Files

```
ui.R        - page layout / inputs / outputs
server.R    - reactive logic: builds the dataset(s), runs distrEstimate(),
              picks the top fits, renders the plots / table / text
global.R    - sources R/main.r (which sources the rest), plus Shiny-only
              helpers: default parameters for the "generate data" panel
              (read from R/distribution.names.txt), a density/pmf lookup,
              dynamic histogram bin selection, and the plotting functions
R/          - your original source files, plus the fixes noted below
```

## Issues: 
1. when hit Generate Numbers, don't overwrite previous numbers. it should generate a new set of numbers. Only when "Reset" button was clicked, all generated data will be deleted. 
2. when new data are generated, it should show up in a new row panel under main panel.
3. for each new set of data, it should have 3 column panes: plot, basic stats, and fit details of distributions that pass the cutoff.
4. In both density plot and fitting details, only show top distribtions that passes the p cutoff, up to 3. 
5. if no distribution pass the cutoff, show the top 3
6. files with multi columns does not work. it should support tab separated .txt, tsv, and csv files. 
7. need to add a "data preview" botton. where when clicked, show data in a separate column panel within the row panel for that set of data. 
8. the histogram frequency and density should use different scale. after fitting, redraw the plot so that the density line won't go out of margin. 

## Issues 2:
1. add a Goodness of fit button below "fitting options" section.
2. when clicked on Generate numbers, or uploaded a file, only show the histogram and basic statistics.
3. when clicked on "Goodness of fit", it will generate the 3 best fitted distribution that pass the p threshold. When no distribution pass the cutoff, show top 3. The density lines should be on top of the histogram. histogram color change from light blue to gray. 
4. show basic statistics as a horizontal box plot (below or above histogram, as you see fit). should label mean and median differently. also show each dots as background, set in light color and be transparent or they don't look dominant. show SD in a way you see approporiate; show min and max data point in different color, and label the values.  
5. add options to copy generated data, and to save data into a file. this option should apply to all generated data. 
6. for data preview. show it in a scrollable text box (not table). No need to show index, just plain text, each data point for a line. 

## issue 3:
1. no need to show the table for basic statistics
2. put the box plot on top of the histogram plot. the two should share the x-axis. 
3. before fitting, the histogram should have a different color (sky blue). 
4. for the box plot, there are overlapps among labelings. try to avoid that. 
5. for the box plot, mean and SD can use the same color. 
6. show data preview on the right side of fitting details. 
7. for copy data/save data, one button should apply to all generated data. 
8. for fitting details, when pvalue is less than 0.0001, show as scientific notation. 
9. for fitting details, put each parameter on its own line. align them properly. 

## Issue 4:
1. copy does not work (doesn't go to clipboard)
2. for copy and save as files, format both as tab separated tables. consider the situation where data have different length. 

# Issue 5
2. add options to remove individual set of data
3. add options of saving fitting results (plot and fitting details) as html. 

## Issue 6
1. seems to refit every dataset when generate new data
2. very slow when fitting discrete (binomial and negative binomial?); expecially when data is large 
3. when fit to bern, estimate prob.
4. p-values are always small with fit with discrete model
5. sometimes it doesn't produce density plot after fitting


## What changed in this round

1. **P-value cutoff reverted to the standard convention.** Default is
   `0.2`; a fit is "good" when `p-value > cutoff`, and the top 3 are the
   ones with the *largest* p-values. (Your instinct was right the first
   time — sorry for the confusion last round.)

2. **One histogram, not two.** There's no separate "preview" plot and
   "fit" plot anymore. Each series gets exactly one histogram. Before
   you run a fit it just shows the plain histogram; after you run a fit
   on that series, the fitted density curves are drawn directly on top
   of that same histogram (same plot, in place) — nothing new is
   created.

3. **Stats beside the histogram, one row per series.**
   - Each numeric series (the one generated vector, or — when you
     upload a file — *every* numeric column in it) gets its own row:
     histogram on the left, a small stats table (n, mean, sd, median,
     min, max) on the right.
   - Uploading a file with multiple numeric columns now previews all of
     them at once, each in its own row, instead of forcing you to pick
     one column just to see it. You still pick exactly one of them
     ("Column / series to fit" in the sidebar) to actually run the
     goodness-of-fit test on.
   - **Save to File** / **Copy to Clipboard** always act on whichever
     series is currently selected to fit, one number per line.

4. **Fixed the "`object 'size' not found`" error for `xnorm`/`xpois`.**
   I couldn't reproduce this in a plain R session — `xnorm`/`xpois`
   don't reference `size` anywhere, and the code that reads generate
   parameters only ever looks up the current distribution's own
   parameter names. My best working theory is a Shiny quirk: **Shiny
   does not clear an input's server-side value just because its UI
   element was removed** (e.g. switching from *binomial*, which has a
   `size` box, to *normal*, which doesn't) — that stale value can
   linger in `input` under its old ID. I couldn't fully confirm this is
   the mechanism without a live browser session, but I closed off the
   whole class of bug regardless: parameter input IDs are now
   **namespaced by distribution** — `genparam_binom_size` vs.
   `genparam_norm_mean` vs. `genparam_pois_lambda`, etc. — so an old
   value under one distribution's ID can now never be read while a
   different distribution is selected, structurally, regardless of the
   exact mechanism. Please let me know if this resolves it or if you
   still see it (and if so, ideally the exact click sequence that
   triggers it, since I can't run the actual browser session here).

## Still true from last round (see previous notes for details)

- `gofTest.chisqC()` renames a fitted weibull's parameters to
  `(shape, scale)` internally — the plot code uses `params$scale`
  directly, not `$lambda`.
- `xgammaFit()` no longer discards your data and fits against a fixed
  synthetic sample.
- Each individual distribution's fit attempt inside `distrEstimate()` is
  wrapped in `tryCatch`, so one distribution failing to optimize (e.g.
  gamma choking on count data containing zeros) doesn't take down the
  whole batch.
- `R/distribution.names.txt` has `default.1`/`default.2`/`default.3`
  columns; `DEFAULT_GEN_PARAMS` in `global.R` is built from those
  instead of being hard-coded.
- Histogram bins scale with sample size (20–50 bins, continuous data).
- Plot heights are in viewport units (`vh`) so they resize with the
  browser window, with extra top margin so titles don't get clipped.

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
