# anovapowersim shinylive app

This repository contains a Shiny app for running `anovapowersim::power_n()` in the browser with shinylive.

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
shinylive::export("app", "docs")
```

Then configure GitHub Pages to serve from the `docs/` folder.
