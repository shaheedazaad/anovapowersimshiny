# anovapowersim shinylive app

This repository contains a Shiny app for running calculation-only ANOVA power analyses in the browser with shinylive. The app vendors the `power_n_calc()`, `power_achieved_calc()`, and `power_sensitivity_calc()` backends, displays the equivalent `anovapowersim` call for reproducible R workflows, and provides PDF reports for every analysis type. (https://shaheedazaad.github.io/anovapowersim/)

## Local Shiny run

```r
old_warn <- getOption("warn")
options(warn = -1)
suppressPackageStartupMessages(library(shiny))
run_app <- shiny::runApp
options(warn = old_warn)
run_app("app")
```

## Build GitHub Pages output

```r
install.packages("shinylive")
shinylive::export(
  "app",
  "docs",
  template_params = list(title = "anovapowersim Shiny app")
)
```

Then configure GitHub Pages to serve from the `docs/` folder.
