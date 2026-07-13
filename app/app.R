suppressWarnings(suppressPackageStartupMessages(library(shiny)))

vendor_candidates <- file.path(
  c(".", dirname(normalizePath("app/app.R", mustWork = FALSE))),
  "vendor",
  "aps",
  "R"
)
vendor_dir <- vendor_candidates[dir.exists(vendor_candidates)][[1]]
for (file in c("utils.R", "power_calc_helpers.R", "power_n_calc.R", "print_methods.R")) {
  source(file.path(vendor_dir, file), local = TRUE)
}
power_n_backend <- power_n_calc

power_n_calc_arg_names <- c(
  "between",
  "within",
  "term",
  "target_pes",
  "power",
  "alpha",
  "n_start",
  "n_max",
  "gpower"
)

power_n_calc_args <- function(args) {
  args[power_n_calc_arg_names]
}

power_n_defaults <- list(
  power = 0.90,
  alpha = 0.05,
  n_start = NULL,
  n_max = 5000L,
  gpower = FALSE
)

power_n_call_args <- function(args) {
  out <- args[c("between", "within", "term", "target_pes")]
  if (!isTRUE(all.equal(args$power, power_n_defaults$power))) {
    out$power <- args$power
  }
  if (!isTRUE(all.equal(args$alpha, power_n_defaults$alpha))) {
    out$alpha <- args$alpha
  }
  if (!is.null(args$n_start)) {
    out$n_start <- args$n_start
  }
  if (!isTRUE(all.equal(args$n_max, power_n_defaults$n_max))) {
    out$n_max <- args$n_max
  }
  if (isTRUE(args$gpower)) {
    out$gpower <- TRUE
  }
  out
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

term_choices <- function(terms) {
  labels <- ifelse(
    grepl(":", terms, fixed = TRUE),
    paste0(terms, " (interaction)"),
    paste0(terms, " (main effect)")
  )
  names(terms) <- labels
  terms
}

term_kind <- function(term) {
  if (grepl(":", term, fixed = TRUE)) "interaction" else "main effect"
}

format_count <- function(x) {
  if (is.na(x)) return("not reached")
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_decimal <- function(x, digits = 3) {
  sub("^0", "", sprintf(paste0("%.", digits, "f"), x))
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
      .progress-panel {
        background: #ecfdf5;
        border: 1px solid #99f6e4;
        border-radius: 8px;
        color: #134e4a;
        margin-bottom: 14px;
        padding: 12px;
      }
      .progress-label {
        font-weight: 650;
        margin-bottom: 8px;
      }
      .progress-detail {
        color: #0f766e;
        font-size: 13px;
        margin: 8px 0 0;
      }
      .progress-log {
        background: #12302d;
        border-radius: 6px;
        color: #d1fae5;
        font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
        font-size: 12px;
        margin: 10px 0 0;
        max-height: 150px;
        overflow: auto;
        padding: 8px;
        white-space: pre-wrap;
      }
      .result-statement {
        background: #eef6f4;
        border: 1px solid #b7d9d3;
        border-left: 6px solid var(--accent);
        border-radius: 8px;
        margin-bottom: 14px;
        padding: 16px 18px;
      }
      .result-sentence {
        font-size: 20px;
        line-height: 1.35;
        margin: 0;
      }
      .result-sentence strong {
        color: var(--accent-dark);
        font-weight: 700;
      }
      .result-statement.not-reached {
        background: #fff1f2;
        border-color: #fecdd3;
        border-left-color: #be123c;
      }
      .result-statement.not-reached .result-sentence strong {
        color: #9f1239;
      }
      .result-detail {
        color: var(--muted);
        margin: 8px 0 0;
      }
      .empty-results {
        background: #fbfcfd;
        border: 1px dashed #b8c2cc;
        border-radius: 8px;
        color: var(--ink);
        padding: 18px;
      }
      .empty-results h3 {
        font-size: 16px;
        line-height: 1.3;
        margin: 0 0 6px;
      }
      .empty-results p {
        color: var(--muted);
        margin: 0;
      }
      .result-view .nav-tabs {
        border-bottom: 1px solid var(--line);
        margin-bottom: 12px;
      }
      .result-view .nav-tabs > li > a {
        color: var(--muted);
        border-radius: 6px 6px 0 0;
        padding: 8px 12px;
      }
      .result-view .nav-tabs > li.active > a,
      .result-view .nav-tabs > li.active > a:focus,
      .result-view .nav-tabs > li.active > a:hover {
        color: var(--accent-dark);
        font-weight: 650;
      }
      .table-scroll {
        overflow-x: auto;
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
        p(class = "subtitle", "Balanced factorial ANOVA power calculation")
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
        " R package. Browser runs use analytic power calculations; the R call below shows the equivalent ",
        code("power_n()"),
        " simulation call for reproducible workflows."
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
      ),
      p(
        "This app is in beta. Please report issues at ",
        tags$a(
          href = "https://github.com/shaheedazaad/anovapowersimshiny",
          target = "_blank",
          rel = "noopener noreferrer",
          "https://github.com/shaheedazaad/anovapowersimshiny"
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
          h2("Power Calculation"),
          selectInput("term", "Term", choices = character()),
          div(
            class = "inline-grid",
            numericInput("target_pes", "Target partial eta squared", value = 0.03, min = 0.001, max = 0.999, step = 0.001),
            numericInput("power", "Target power", value = 0.90, min = 0.001, max = 0.999, step = 0.01),
            numericInput("alpha", "Alpha", value = 0.05, min = 0.001, max = 0.999, step = 0.001),
            numericInput("n_max", "Search up to n per cell", value = power_n_defaults$n_max, min = 1, step = 100)
          ),
          tags$details(
            class = "advanced-settings",
            tags$summary("Advanced settings"),
            div(
              class = "checks",
              checkboxInput("use_n_start", "Set n_start", value = FALSE),
              conditionalPanel(
                condition = "input.use_n_start",
                numericInput("n_start", "n_start", value = 10, min = 1, step = 1)
              ),
              checkboxInput("gpower", "G*Power convention", value = FALSE)
            )
          ),
          div(
            class = "run-row",
            actionButton("run", "Calculate", icon = icon("play"), class = "btn-primary"),
            uiOutput("run_status")
          )
        )
      ),
      div(
        div(
          class = "panel",
          h2("Results"),
          uiOutput("progress"),
          uiOutput("summary"),
          uiOutput("warnings"),
          uiOutput("result_view_ui")
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
  progress_state <- reactiveVal(list(
    active = FALSE,
    lines = character()
  ))

  set_progress_state <- function(state) {
    progress_state(state)
  }

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
      term_choices(design_terms(df$name))
    } else {
      character()
    }
    values <- unname(choices)
    selected <- isolate(input$term)
    if (!length(selected) || !selected %in% values) {
      selected <- if (length(values)) values[[1]] else character()
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
      alpha = input$alpha,
      n_start = if (isTRUE(input$use_n_start)) as.integer(input$n_start) else NULL,
      n_max = as.integer(input$n_max),
      gpower = isTRUE(input$gpower)
    )
    args
  })

  output$call <- renderText({
    tryCatch(format_call(power_n_call_args(power_args())), error = function(e) "")
  })

  observeEvent(input$run, {
    args <- power_n_calc_args(power_args())
    captured_warnings <- character()
    run_warnings(character())
    result(NULL)
    running(TRUE)
    set_progress_state(list(
      active = TRUE,
      lines = character()
    ))

    session$onFlushed(function() {
      on.exit({
        running(FALSE)
        set_progress_state(list(
          active = FALSE,
          lines = character()
        ))
      }, add = TRUE)

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
            captured_warnings <<- unique(c(captured_warnings, msg))
            run_warnings(captured_warnings)
          }
          restart <- findRestart("muffleWarning")
          if (!is.null(restart)) {
            invokeRestart(restart)
          }
        }
      )
      result(fit)
    }, once = TRUE)
  })

  output$run_status <- renderUI({
    if (isTRUE(running())) {
      span(class = "text-muted", "Calculating...")
    } else {
      NULL
    }
  })

  output$progress <- renderUI({
    state <- progress_state()
    if (!isTRUE(state$active) && !length(state$lines)) {
      return(NULL)
    }
    div(
      class = "progress-panel",
      div(
        class = "progress-label",
        "Power calculation is running"
      ),
      div(
        class = "progress-detail",
        "Keep this browser tab open until the results appear."
      )
    )
  })

  output$summary <- renderUI({
    fit <- result()
    if (is.null(fit)) {
      return(div(
        class = "empty-results",
        h3("Ready to calculate power"),
        p("Enter the design details, choose the ANOVA term and target effect size, then click Calculate.")
      ))
    }
    kind <- term_kind(fit$term)
    power_text <- format_decimal(fit$power)
    pes_text <- format_decimal(fit$target_pes)
    alpha_text <- format_decimal(fit$alpha)
    n_per_cell_text <- format_count(fit$n_needed)
    total_n_text <- format_count(fit$total_n_needed)
    between_cells <- fit$design$n_between_cells %||% 1L
    detail <- if (is.na(fit$n_needed)) {
      tagList(
        sprintf(
          "The search ran up to n_max = %s participants per between-subjects cell at alpha = %s with target partial eta squared = %s.",
          format_count(max(fit$results$n_per_cell, na.rm = TRUE)),
          alpha_text,
          pes_text
        ),
        " Increase ",
        strong("Search up to n per cell"),
        " and calculate again to extend the search."
      )
    } else if (between_cells > 1L) {
      sprintf(
        "This corresponds to %s participants per between-subjects cell at alpha = %s with target partial eta squared = %s.",
        n_per_cell_text,
        alpha_text,
        pes_text
      )
    } else {
      sprintf(
        "This is based on alpha = %s and target partial eta squared = %s.",
        alpha_text,
        pes_text
      )
    }

    div(
      class = if (is.na(fit$n_needed)) "result-statement not-reached" else "result-statement",
      if (is.na(fit$n_needed)) {
        p(
          class = "result-sentence",
          "The target power of ",
          strong(power_text),
          " was not reached for the ",
          strong(fit$term),
          " ",
          kind,
          "."
        )
      } else {
        p(
          class = "result-sentence",
          strong(total_n_text),
          " participants are needed to achieve ",
          strong(power_text),
          " power for the ",
          strong(fit$term),
          " ",
          kind,
          "."
        )
      },
      p(class = "result-detail", detail)
    )
  })

  output$result_view_ui <- renderUI({
    req(result())
    div(
      class = "result-view",
      tabsetPanel(
        id = "result_view",
        type = "tabs",
        tabPanel("Plot", plotOutput("curve_plot", height = 300)),
        tabPanel("Table", div(class = "table-scroll", tableOutput("results_table")))
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
    validate(need(nrow(results) > 0, "No power rows returned."))
    validate(need(any(is.finite(results$power_calc)), "No finite calculated power values returned."))

    has_sim <- "power_sim" %in% names(results) && any(is.finite(results$power_sim))
    ylim <- range(c(results$power_calc, if (has_sim) results$power_sim, fit$power), na.rm = TRUE)
    ylim <- c(max(0, ylim[[1]] - 0.05), min(1, ylim[[2]] + 0.05))
    plot(
      results$n_per_cell,
      results$power_calc,
      type = "b",
      pch = 17,
      col = "#0f766e",
      ylim = ylim,
      xlab = "n per between-subjects cell",
      ylab = "Power"
    )
    if (has_sim) {
      lines(results$n_per_cell, results$power_sim, type = "b", pch = 19, col = "#374151")
    }
    abline(h = fit$power, lty = 2, col = "#9f1239")
    legend_labels <- c("Calculated", if (has_sim) "Simulated", "Target")
    legend_cols <- c("#0f766e", if (has_sim) "#374151", "#9f1239")
    legend_pch <- c(17, if (has_sim) 19, NA)
    legend_lty <- c(1, if (has_sim) 1, 2)
    legend("bottomright", legend = legend_labels, col = legend_cols,
           lty = legend_lty, pch = legend_pch, bty = "n")
  })

  output$results_table <- renderTable({
    fit <- result()
    req(fit)
    results <- as.data.frame(fit$results)
    results[, setdiff(names(results), c("n_sims", "power_sim")), drop = FALSE]
  }, striped = TRUE, bordered = TRUE, digits = 3)
}

shinyApp(ui, server)
