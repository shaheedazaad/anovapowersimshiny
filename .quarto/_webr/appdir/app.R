suppressWarnings(suppressPackageStartupMessages(library(shiny)))

anovapowersim_pkg <- paste0("anova", "powersim")
if (suppressWarnings(requireNamespace(anovapowersim_pkg, quietly = TRUE))) {
  power_n_backend <- getExportedValue(anovapowersim_pkg, "power_n")
} else {
  vendor_candidates <- file.path(
    c(".", dirname(normalizePath("app/app.R", mustWork = FALSE))),
    "vendor",
    "anovapowersim",
    "R"
  )
  vendor_dir <- vendor_candidates[dir.exists(vendor_candidates)][[1]]
  for (file in c("utils.R", "compute_scale_factor.R", "power_curve.R", "print_methods.R")) {
    source(file.path(vendor_dir, file), local = TRUE)
  }
  power_n_backend <- power_n
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

should_show_power_warning <- function(message) {
  hidden_patterns <- c(
    "^package .+ was built under R version",
    "^Package version mismatch for ",
    "^namespace .+ is not available and has been replaced",
    "^replacing previous import"
  )
  !any(grepl(paste(hidden_patterns, collapse = "|"), message))
}

factor_row <- function(id, type, default_name, default_levels) {
  div(
    class = "factor-row",
    selectInput(
      inputId = paste0("factor_type_", id),
      label = NULL,
      choices = c("Between" = "between", "Within" = "within"),
      selected = type,
      width = "100%"
    ),
    textInput(
      inputId = paste0("factor_name_", id),
      label = NULL,
      value = default_name,
      placeholder = "name",
      width = "100%"
    ),
    numericInput(
      inputId = paste0("factor_levels_", id),
      label = NULL,
      value = default_levels,
      min = 2,
      step = 1,
      width = "100%"
    ),
    actionButton(
      inputId = paste0("remove_factor_", id),
      label = "",
      icon = icon("trash"),
      class = "btn-outline-danger remove-factor",
      onclick = sprintf("Shiny.setInputValue('remove_factor', %s, {priority: 'event'})", id)
    )
  )
}

clean_factor_name <- function(x) {
  x <- trimws(x %||% "")
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x) || grepl("^[0-9]", x)) {
    return(NA_character_)
  }
  x
}

named_counts <- function(df, type) {
  rows <- df[df$type == type, , drop = FALSE]
  if (!nrow(rows)) {
    return(NULL)
  }
  out <- as.integer(rows$levels)
  names(out) <- rows$name
  out
}

design_terms <- function(factor_names) {
  if (!length(factor_names)) {
    return(character())
  }
  unlist(
    lapply(seq_along(factor_names), function(k) {
      combn(factor_names, k, FUN = paste, collapse = ":", simplify = TRUE)
    }),
    use.names = FALSE
  )
}

format_call <- function(args) {
  value_to_code <- function(x) {
    if (is.null(x)) {
      return("NULL")
    }
    if (is.logical(x)) {
      return(if (isTRUE(x)) "TRUE" else "FALSE")
    }
    if (is.character(x)) {
      return(sprintf('"%s"', x))
    }
    if (length(x) > 1 || !is.null(names(x))) {
      pieces <- paste(sprintf("%s = %s", names(x), x), collapse = ", ")
      return(sprintf("c(%s)", pieces))
    }
    as.character(x)
  }

  lines <- sprintf("  %s = %s", names(args), vapply(args, value_to_code, character(1)))
  paste0("power_n(\n", paste(lines, collapse = ",\n"), "\n)")
}

ui <- fluidPage(
  tags$head(
    tags$title("anovapowersim Shiny app"),
    tags$style(HTML("
      :root {
        --ink: #1f2937;
        --muted: #6b7280;
        --line: #d7dde5;
        --panel: #ffffff;
        --bg: #f5f7fa;
        --accent: #0f766e;
        --accent-dark: #115e59;
      }
      body {
        background: var(--bg);
        color: var(--ink);
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
      }
      .app-shell {
        max-width: 1280px;
        margin: 0 auto;
        padding: 24px;
      }
      .app-header {
        display: flex;
        justify-content: space-between;
        gap: 20px;
        align-items: end;
        border-bottom: 1px solid var(--line);
        padding-bottom: 18px;
        margin-bottom: 18px;
      }
      h1 {
        font-size: 28px;
        line-height: 1.15;
        margin: 0 0 6px;
        letter-spacing: 0;
      }
      .subtitle {
        color: var(--muted);
        margin: 0;
      }
      .layout {
        display: grid;
        grid-template-columns: minmax(330px, 430px) minmax(0, 1fr);
        gap: 18px;
        align-items: start;
      }
      .panel {
        background: var(--panel);
        border: 1px solid var(--line);
        border-radius: 8px;
        padding: 16px;
      }
      .panel + .panel {
        margin-top: 14px;
      }
      .panel h2 {
        font-size: 16px;
        margin: 0 0 12px;
        letter-spacing: 0;
      }
      .info-panel {
        background: #eef6f4;
        border: 1px solid #b7d9d3;
        border-radius: 8px;
        margin-bottom: 18px;
        padding: 14px 16px;
      }
      .info-panel p {
        margin: 0;
      }
      .info-panel p + p {
        margin-top: 8px;
      }
      .info-panel a {
        color: var(--accent-dark);
        font-weight: 600;
      }
      .factor-row {
        display: grid;
        grid-template-columns: minmax(92px, .85fr) minmax(96px, 1fr) minmax(70px, .65fr) 38px;
        gap: 8px;
        align-items: start;
        margin-bottom: 8px;
      }
      .factor-header {
        color: var(--muted);
        display: grid;
        font-size: 12px;
        font-weight: 650;
        gap: 8px;
        grid-template-columns: minmax(92px, .85fr) minmax(96px, 1fr) minmax(70px, .65fr) 38px;
        margin: 0 0 6px;
        text-transform: uppercase;
      }
      .factor-row .form-group,
      .inline-grid .form-group {
        margin-bottom: 0;
      }
      .factor-row .form-control,
      .factor-row .selectize-control,
      .factor-row .selectize-input {
        min-width: 0;
        width: 100%;
      }
      .remove-factor {
        align-items: center;
        display: inline-flex;
        height: 38px;
        justify-content: center;
        padding: 0;
        width: 38px;
      }
      .inline-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 12px;
      }
      .advanced-settings {
        border-top: 1px solid var(--line);
        margin-top: 14px;
        padding-top: 12px;
      }
      .advanced-settings summary {
        color: var(--accent-dark);
        cursor: pointer;
        font-weight: 600;
        margin-bottom: 10px;
      }
      .advanced-settings[open] summary {
        margin-bottom: 12px;
      }
      .checks {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 0 12px;
      }
      .run-row {
        display: flex;
        gap: 10px;
        align-items: center;
        margin-top: 14px;
      }
      .btn-primary,
      .btn-primary:focus {
        background: var(--accent);
        border-color: var(--accent);
      }
      .btn-primary:hover {
        background: var(--accent-dark);
        border-color: var(--accent-dark);
      }
      .summary-grid {
        display: grid;
        grid-template-columns: repeat(3, minmax(0, 1fr));
        gap: 12px;
        margin-bottom: 14px;
      }
      .metric {
        border: 1px solid var(--line);
        border-radius: 8px;
        padding: 12px;
        background: #fbfcfd;
      }
      .metric-label {
        color: var(--muted);
        font-size: 12px;
        text-transform: uppercase;
        letter-spacing: .04em;
      }
      .metric-value {
        font-size: 24px;
        font-weight: 650;
        margin-top: 4px;
      }
      .code-box {
        background: #111827;
        color: #e5e7eb;
        border-radius: 8px;
        padding: 12px;
        overflow: auto;
        font-size: 13px;
      }
      .table {
        font-size: 13px;
      }
      .validation {
        color: #9f1239;
        margin: 10px 0 0;
      }
      .warnings-panel {
        background: #fff7ed;
        border: 1px solid #fed7aa;
        border-radius: 8px;
        color: #7c2d12;
        margin-bottom: 14px;
        padding: 12px;
      }
      .warnings-panel h3 {
        font-size: 14px;
        margin: 0 0 8px;
      }
      .warnings-panel ul {
        margin: 0;
        padding-left: 18px;
      }
      .warnings-panel li + li {
        margin-top: 6px;
      }
      .shiny-output-error-validation {
        color: #9f1239;
      }
      @media (max-width: 900px) {
        .app-shell { padding: 16px; }
        .app-header { display: block; }
        .layout { grid-template-columns: 1fr; }
        .summary-grid { grid-template-columns: 1fr; }
      }
      @media (max-width: 560px) {
        .factor-row {
          grid-template-columns: minmax(0, 1fr) minmax(0, 1fr) 70px 38px;
        }
        .factor-header {
          grid-template-columns: minmax(0, 1fr) minmax(0, 1fr) 70px 38px;
        }
        .inline-grid,
        .checks {
          grid-template-columns: 1fr;
        }
      }
    "))
  ),
  div(
    class = "app-shell",
    div(
      class = "app-header",
      div(
        h1("anovapowersim Shiny app"),
        p(class = "subtitle", "Balanced factorial ANOVA sample size search")
      ),
      tags$a(
        href = "https://shaheedazaad.github.io/anovapowersim/reference/power_n.html",
        target = "_blank",
        rel = "noopener noreferrer",
        "power_n() reference"
      )
    ),
    div(
      class = "info-panel",
      p(
        "This app is a browser interface to the ",
        tags$a(
          href = "https://shaheedazaad.github.io/anovapowersim/",
          target = "_blank",
          rel = "noopener noreferrer",
          "anovapowersim"
        ),
        " R package. For the full feature set and reproducible workflows, use the R package directly."
      ),
      p(
        "Citation: Azaad, S. (2026). A priori power analysis for ANOVA interaction effects with the anovapowersim R package: a short introduction. ",
        tags$a(
          href = "https://doi.org/10.31234/osf.io/86rsy_v1",
          target = "_blank",
          rel = "noopener noreferrer",
          "https://doi.org/10.31234/osf.io/86rsy_v1"
        ),
        "."
      )
    ),
    div(
      class = "layout",
      div(
        div(
          class = "panel",
          h2("Design"),
          div(
            class = "factor-header",
            span("Type"),
            span("Factor"),
            span("Levels"),
            span("")
          ),
          uiOutput("factor_rows"),
          actionButton("add_factor", "Add factor", icon = icon("plus"), class = "btn-outline-secondary"),
          uiOutput("factor_validation")
        ),
        div(
          class = "panel",
          h2("Power Analysis"),
          selectInput("term", "Term", choices = character()),
          div(
            class = "inline-grid",
            numericInput("target_pes", "Target partial eta squared", value = 0.03, min = 0.001, max = 0.999, step = 0.001),
            numericInput("power", "Target power", value = 0.90, min = 0.001, max = 0.999, step = 0.01),
            numericInput("n_sims", "Simulations per sample size", value = 10000, min = 1, step = 100),
            numericInput("alpha", "Alpha", value = 0.05, min = 0.001, max = 0.999, step = 0.001),
            numericInput("seed", "Seed", value = 123, min = 1, step = 1)
          ),
          div(
            class = "checks",
            checkboxInput("use_seed", "Set seed", value = TRUE),
            checkboxInput("parallel", "Parallel", value = TRUE)
          ),
          tags$details(
            class = "advanced-settings",
            tags$summary("Advanced settings"),
            div(
              class = "inline-grid",
              selectInput("ss_type", "Sums of squares", choices = c("I", "II", "III"), selected = "III"),
              numericInput("n_max", "Maximum n per cell", value = 1000, min = 1, step = 1),
              numericInput("tol", "Tolerance", value = 0.03, min = 0.001, step = 0.001)
            ),
            div(
              class = "checks",
              checkboxInput("use_n_start", "Set n_start", value = FALSE),
              conditionalPanel(
                condition = "input.use_n_start",
                numericInput("n_start", "n_start", value = 10, min = 1, step = 1)
              ),
              checkboxInput("use_cores", "Set cores", value = FALSE),
              conditionalPanel(
                condition = "input.use_cores",
                numericInput("cores", "cores", value = 2, min = 1, step = 1)
              ),
              checkboxInput("gpower", "G*Power convention", value = FALSE),
              checkboxInput("progress", "Text progress", value = FALSE)
            )
          ),
          div(
            class = "run-row",
            actionButton("run", "Run", icon = icon("play"), class = "btn-primary"),
            uiOutput("run_status")
          )
        )
      ),
      div(
        div(
          class = "panel",
          h2("Results"),
          uiOutput("summary"),
          uiOutput("warnings"),
          plotOutput("curve_plot", height = 300),
          tableOutput("results_table")
        ),
        div(
          class = "panel",
          h2("R Call"),
          verbatimTextOutput("call", placeholder = TRUE)
        )
      )
    )
  )
)

server <- function(input, output, session) {
  factors <- reactiveVal(data.frame(
    id = 1:2,
    type = c("between", "within"),
    default_name = c("group", "time"),
    default_levels = c(2L, 3L),
    stringsAsFactors = FALSE
  ))
  next_factor_id <- reactiveVal(3L)
  running <- reactiveVal(FALSE)
  result <- reactiveVal(NULL)
  run_warnings <- reactiveVal(character())

  observeEvent(input$add_factor, {
    id <- next_factor_id()
    current <- factors()
    factors(rbind(
      current,
      data.frame(
        id = id,
        type = "between",
        default_name = paste0("factor", id),
        default_levels = 2L,
        stringsAsFactors = FALSE
      )
    ))
    next_factor_id(id + 1L)
  })

  observeEvent(input$remove_factor, {
    id <- input$remove_factor
    factors(factors()[factors()$id != id, , drop = FALSE])
  })

  current_design <- reactive({
    rows <- factors()
    if (!nrow(rows)) {
      return(data.frame(id = integer(), type = character(), name = character(), levels = integer()))
    }
    rows$type <- vapply(rows$id, function(id) {
      input[[paste0("factor_type_", id)]] %||% rows$type[rows$id == id]
    }, character(1))
    rows$name <- vapply(rows$id, function(id) {
      clean_factor_name(input[[paste0("factor_name_", id)]] %||% rows$default_name[rows$id == id])
    }, character(1))
    rows$levels <- vapply(rows$id, function(id) {
      value <- input[[paste0("factor_levels_", id)]] %||% rows$default_levels[rows$id == id]
      as.integer(value)
    }, integer(1))
    rows[, c("id", "type", "name", "levels"), drop = FALSE]
  })

  design_problem <- reactive({
    df <- current_design()
    if (!nrow(df)) {
      return("Add at least one factor.")
    }
    if (any(is.na(df$name))) {
      return("Factor names must start with a letter or underscore and use letters, numbers, or underscores.")
    }
    if (anyDuplicated(df$name)) {
      return("Factor names must be unique.")
    }
    if (any(is.na(df$levels) | df$levels < 2L)) {
      return("Each factor needs at least two levels.")
    }
    NULL
  })

  output$factor_rows <- renderUI({
    rows <- factors()
    if (!nrow(rows)) {
      return(NULL)
    }
    tagList(lapply(seq_len(nrow(rows)), function(i) {
      factor_row(
        id = rows$id[[i]],
        type = rows$type[[i]],
        default_name = rows$default_name[[i]],
        default_levels = rows$default_levels[[i]]
      )
    }))
  })

  output$factor_validation <- renderUI({
    problem <- design_problem()
    if (is.null(problem)) {
      return(NULL)
    }
    div(class = "validation", problem)
  })

  observe({
    problem <- design_problem()
    df <- current_design()
    choices <- if (is.null(problem)) {
      design_terms(df$name)
    } else {
      character()
    }
    selected <- isolate(input$term)
    if (!selected %in% choices) {
      selected <- if (length(choices)) choices[[1]] else character()
    }
    updateSelectInput(session, "term", choices = choices, selected = selected)
  })

  power_args <- reactive({
    problem <- design_problem()
    validate(need(is.null(problem), problem))

    df <- current_design()
    req(input$term)
    args <- list(
      between = named_counts(df, "between"),
      within = named_counts(df, "within"),
      term = input$term,
      target_pes = input$target_pes,
      power = input$power,
      n_sims = as.integer(input$n_sims),
      alpha = input$alpha,
      ss_type = input$ss_type,
      n_start = if (isTRUE(input$use_n_start)) as.integer(input$n_start) else NULL,
      n_max = as.integer(input$n_max),
      tol = input$tol,
      gpower = isTRUE(input$gpower),
      progress = isTRUE(input$progress),
      parallel = isTRUE(input$parallel),
      cores = if (isTRUE(input$use_cores)) as.integer(input$cores) else NULL,
      seed = if (isTRUE(input$use_seed)) as.integer(input$seed) else NULL
    )
    args
  })

  output$call <- renderText({
    tryCatch(format_call(power_args()), error = function(e) "")
  })

  observeEvent(input$run, {
    args <- power_args()
    run_warnings(character())
    running(TRUE)
    on.exit(running(FALSE), add = TRUE)

    withProgress(message = "Running simulations", value = 0.2, {
      fit <- withCallingHandlers(
        tryCatch(
          do.call(power_n_backend, args),
          error = function(e) {
            showNotification(conditionMessage(e), type = "error", duration = 10)
            NULL
          }
        ),
        warning = function(w) {
          msg <- conditionMessage(w)
          if (should_show_power_warning(msg)) {
            run_warnings(unique(c(run_warnings(), msg)))
          }
          invokeRestart("muffleWarning")
        }
      )
      incProgress(0.8)
      result(fit)
    })
  })

  output$run_status <- renderUI({
    if (isTRUE(running())) {
      span(class = "text-muted", "Running...")
    } else {
      NULL
    }
  })

  output$summary <- renderUI({
    fit <- result()
    if (is.null(fit)) {
      return(div(class = "text-muted", "No run yet."))
    }
    div(
      class = "summary-grid",
      div(
        class = "metric",
        div(class = "metric-label", "n per cell"),
        div(class = "metric-value", ifelse(is.na(fit$n_needed), "NA", fit$n_needed))
      ),
      div(
        class = "metric",
        div(class = "metric-label", "total N"),
        div(class = "metric-value", ifelse(is.na(fit$total_n_needed), "NA", fit$total_n_needed))
      ),
      div(
        class = "metric",
        div(class = "metric-label", "term"),
        div(class = "metric-value", fit$term)
      )
    )
  })

  output$warnings <- renderUI({
    warnings <- run_warnings()
    if (!length(warnings)) {
      return(NULL)
    }
    div(
      class = "warnings-panel",
      h3("Warnings"),
      tags$ul(lapply(warnings, tags$li))
    )
  })

  output$curve_plot <- renderPlot({
    fit <- result()
    req(fit)
    results <- as.data.frame(fit$results)
    validate(need(nrow(results) > 0, "No simulation rows returned."))

    ylim <- range(c(results$power_calc, results$power_sim, fit$power), na.rm = TRUE)
    ylim <- c(max(0, ylim[[1]] - 0.05), min(1, ylim[[2]] + 0.05))
    plot(
      results$n_per_cell,
      results$power_sim,
      type = "b",
      pch = 19,
      col = "#0f766e",
      ylim = ylim,
      xlab = "n per between-subjects cell",
      ylab = "Power"
    )
    lines(results$n_per_cell, results$power_calc, type = "b", pch = 17, col = "#374151")
    abline(h = fit$power, lty = 2, col = "#9f1239")
    legend(
      "bottomright",
      legend = c("Simulated", "Calculated", "Target"),
      col = c("#0f766e", "#374151", "#9f1239"),
      lty = c(1, 1, 2),
      pch = c(19, 17, NA),
      bty = "n"
    )
  })

  output$results_table <- renderTable({
    fit <- result()
    req(fit)
    as.data.frame(fit$results)
  }, striped = TRUE, bordered = TRUE, digits = 3)
}

shinyApp(ui, server)

