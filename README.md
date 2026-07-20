# anovapowersim shinylive app

This repository contains a Shiny app for calculating ANOVA power in the browser with shinylive. The app uses the vendored calculation-only `power_n_calc()` implementation and displays an equivalent `anovapowersim::power_n_calc()` call for reproducible R workflows. (https://shaheedazaad.github.io/anovapowersim/)

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
