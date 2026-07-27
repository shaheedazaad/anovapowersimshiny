suppressWarnings(suppressPackageStartupMessages(library(shiny)))

vendor_candidates <- file.path(
  c(".", dirname(normalizePath("app/app.R", mustWork = FALSE))),
  "vendor",
  "aps",
  "R"
)
vendor_dir <- vendor_candidates[dir.exists(vendor_candidates)][[1]]
for (file in c(
  "utils.R",
  "power_calc_helpers.R",
  "power_n_calc.R",
  "power_effect_calc.R",
  "print_methods.R"
)) {
  source(file.path(vendor_dir, file), local = TRUE)
}
power_n_backend <- power_n_calc
power_achieved_backend <- power_achieved_calc
power_sensitivity_backend <- power_sensitivity_calc

power_n_calc_arg_names <- c(
  "between",
  "within",
  "term",
  "target_pes",
  "power",
  "alpha",
  "n_start",
  "n_max",
  "gpower",
  "epsilon"
)

power_n_calc_args <- function(args) {
  args[power_n_calc_arg_names]
}

power_n_defaults <- list(
  power = 0.90,
  alpha = 0.05,
  n_start = NULL,
  n_max = 5000L,
  gpower = FALSE,
  epsilon = 1
)

power_effect_defaults <- list(
  power = 0.90,
  alpha = 0.05,
  pes_min = 1e-6,
  pes_max = 0.99,
  pes_tol = 0.001,
  gpower = FALSE,
  epsilon = 1
)

power_n_call_args <- function(args) {
  out <- args[c("between", "within", "term", "target_pes", "epsilon")]
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

power_achieved_call_args <- function(args) {
  out <- args[c("between", "within", "term", "target_pes", "n", "epsilon")]
  if (!isTRUE(all.equal(args$alpha, power_effect_defaults$alpha))) {
    out$alpha <- args$alpha
  }
  if (isTRUE(args$gpower)) {
    out$gpower <- TRUE
  }
  out
}

power_sensitivity_call_args <- function(args) {
  out <- args[c("between", "within", "term", "n", "epsilon")]
  for (field in c("power", "alpha", "pes_min", "pes_max", "pes_tol")) {
    if (!isTRUE(all.equal(args[[field]], power_effect_defaults[[field]]))) {
      out[[field]] <- args[[field]]
    }
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

preserve_factor_inputs <- function(rows, input) {
  if (!nrow(rows)) {
    return(rows)
  }
  rows$type <- vapply(rows$id, function(id) {
    input[[paste0("factor_type_", id)]] %||% rows$type[rows$id == id]
  }, character(1))
  rows$default_name <- vapply(rows$id, function(id) {
    input[[paste0("factor_name_", id)]] %||% rows$default_name[rows$id == id]
  }, character(1))
  rows$default_levels <- vapply(rows$id, function(id) {
    value <- input[[paste0("factor_levels_", id)]] %||%
      rows$default_levels[rows$id == id]
    as.integer(value)
  }, integer(1))
  rows
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

format_call <- function(function_name, args) {
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
  paste0(
    "install.packages(\"pak\")\n",
    "pak::pak(\"shaheedazaad/anovapowersim\")\n",
    "library(anovapowersim)\n\n",
    function_name, "(\n", paste(lines, collapse = ",\n"), "\n)"
  )
}

draw_power_curve <- function(fit, show_x_label = TRUE) {
  results <- as.data.frame(fit$results)
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
    xlab = if (show_x_label) "n per between-subjects cell" else "",
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
}

draw_sensitivity_curve <- function(fit, legend_cex = 1) {
  results <- as.data.frame(fit$results)
  plot(
    results$target_pes,
    results$power_calc,
    type = "b",
    pch = 17,
    col = "#0f766e",
    ylim = c(0, 1),
    xlab = "Partial eta squared",
    ylab = "Power"
  )
  abline(h = fit$power, lty = 2, col = "#9f1239")
  if (is.finite(fit$pes_needed)) {
    abline(v = fit$pes_needed, lty = 3, col = "#374151")
  }
  legend(
    "bottomright",
    legend = c(
      "Calculated",
      "Target power",
      if (is.finite(fit$pes_needed)) "Detectable effect"
    ),
    col = c(
      "#0f766e",
      "#9f1239",
      if (is.finite(fit$pes_needed)) "#374151"
    ),
    lty = c(1, 2, if (is.finite(fit$pes_needed)) 3),
    pch = c(17, NA, if (is.finite(fit$pes_needed)) NA),
    cex = legend_cex,
    bty = "n"
  )
}

report_factor_text <- function(design, type) {
  names <- design[[type]] %||% character()
  level_counts <- design$level_counts %||% integer()
  if (!length(names)) {
    return("None")
  }
  paste(
    sprintf("%s (%s levels)", names, level_counts[names]),
    collapse = "; "
  )
}

report_results_table <- function(fit) {
  results <- as.data.frame(fit$results)
  fields <- c("n_per_cell", "total_n", "num_df", "den_df", "ncp", "power_calc")
  labels <- c("n per cell", "Total N", "Numerator df", "Denominator df", "NCP", "Power")
  table <- results[, fields, drop = FALSE]
  names(table) <- labels
  table[["n per cell"]] <- format(table[["n per cell"]], trim = TRUE, scientific = FALSE)
  table[["Total N"]] <- format(table[["Total N"]], trim = TRUE, scientific = FALSE)
  for (field in c("Numerator df", "Denominator df", "NCP", "Power")) {
    table[[field]] <- formatC(table[[field]], format = "f", digits = 3)
  }
  as.matrix(table)
}

report_sensitivity_table <- function(fit) {
  results <- as.data.frame(fit$results)
  fields <- c("target_pes", "num_df", "den_df", "ncp", "power_calc")
  labels <- c("Partial eta sq.", "Num. df", "Den. df", "NCP", "Power")
  table <- results[, fields, drop = FALSE]
  names(table) <- labels
  table[["Partial eta sq."]] <- formatC(
    table[["Partial eta sq."]],
    format = "f",
    digits = 6
  )
  for (field in c("Num. df", "Den. df", "NCP", "Power")) {
    table[[field]] <- formatC(table[[field]], format = "f", digits = 3)
  }
  as.matrix(table)
}

report_fixed_n_text <- function(fit) {
  if (length(fit$design$between)) {
    sprintf(
      "Fixed sample: %s per between-subjects cell (%s total participants)",
      format_count(fit$n),
      format_count(fit$total_n)
    )
  } else {
    sprintf("Fixed sample: %s total participants", format_count(fit$total_n))
  }
}

report_convention_text <- function(fit) {
  sprintf(
    "Nonsphericity correction: %.3f    G*Power convention: %s",
    fit$epsilon,
    if (isTRUE(fit$gpower)) "yes" else "no"
  )
}

draw_report_footer <- function() {
  graphics::text(
    0.5,
    0.02,
    "Generated by the anovapowersim Shiny app",
    cex = 0.6,
    col = "#6b7280"
  )
}

draw_report_table <- function(table, title = "Calculated-power results") {
  x <- seq(0.08, 0.92, length.out = ncol(table))
  plot.new()
  text(0.5, 0.96, title, font = 2, cex = 0.85)
  text(x, 0.88, colnames(table), font = 2, cex = 0.64)
  segments(0.04, 0.85, 0.96, 0.85, col = "#6b7280")
  row_spacing <- min(0.04, 0.70 / max(1L, nrow(table) - 1L))
  y <- 0.79 - (seq_len(nrow(table)) - 1L) * row_spacing
  for (i in seq_len(nrow(table))) {
    text(x, y[[i]], table[i, ], cex = 0.64)
  }
}

write_power_report <- function(file, fit) {
  table <- report_results_table(fit)

  grDevices::pdf(file, width = 8.5, height = 11)
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::text(0.5, 0.96, "ANOVA power analysis report", font = 2, cex = 1.25)
  graphics::text(0.06, 0.91, paste("Between-subjects factors:", report_factor_text(fit$design, "between")), adj = 0, cex = 0.8)
  graphics::text(0.06, 0.875, paste("Within-subjects factors:", report_factor_text(fit$design, "within")), adj = 0, cex = 0.8)
  graphics::text(0.06, 0.84, paste("Tested term:", fit$term), adj = 0, cex = 0.8)
  graphics::text(
    0.06, 0.805,
    sprintf("Target power: %.3f    Alpha: %.3f    Target partial eta squared: %.3f", fit$power, fit$alpha, fit$target_pes),
    adj = 0, cex = 0.8
  )
  graphics::text(0.06, 0.775, sprintf("Nonsphericity correction: %.3f", fit$epsilon), adj = 0, cex = 0.8)
  result_text <- if (is.na(fit$n_needed)) {
    "Target power was not reached within the requested search range."
  } else {
    sprintf("Required sample: %s per between-subjects cell (%s total participants).", fit$n_needed, fit$total_n_needed)
  }
  graphics::text(0.06, 0.745, result_text, adj = 0, font = 2, cex = 0.8)
  draw_report_footer()
  graphics::par(fig = c(0.11, 0.91, 0.48, 0.70), new = TRUE, mar = c(3.5, 3.8, 0.5, 0.5))
  draw_power_curve(fit, show_x_label = FALSE)
  graphics::par(fig = c(0.06, 0.94, 0.06, 0.45), new = TRUE, mar = c(0, 0, 0, 0))
  draw_report_table(table, "Sample-size search results")
}

write_achieved_report <- function(file, fit) {
  table <- report_results_table(fit)

  grDevices::pdf(file, width = 8.5, height = 11)
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::text(0.5, 0.96, "ANOVA achieved-power report", font = 2, cex = 1.25)
  graphics::text(0.06, 0.91, paste("Between-subjects factors:", report_factor_text(fit$design, "between")), adj = 0, cex = 0.8)
  graphics::text(0.06, 0.875, paste("Within-subjects factors:", report_factor_text(fit$design, "within")), adj = 0, cex = 0.8)
  graphics::text(0.06, 0.84, paste("Tested term:", fit$term), adj = 0, cex = 0.8)
  graphics::text(0.06, 0.805, report_fixed_n_text(fit), adj = 0, cex = 0.8)
  graphics::text(
    0.06, 0.77,
    sprintf("Alpha: %.3f    Partial eta squared: %.4f", fit$alpha, fit$target_pes),
    adj = 0, cex = 0.8
  )
  graphics::text(0.06, 0.735, report_convention_text(fit), adj = 0, cex = 0.8)
  graphics::rect(0.06, 0.59, 0.94, 0.68, col = "#b7d9d3", border = "#0f766e")
  graphics::text(
    0.5, 0.635,
    sprintf("Achieved power: %.3f", fit$achieved_power),
    font = 2, cex = 1.15, col = "#115e59"
  )
  draw_report_footer()
  graphics::par(fig = c(0.06, 0.94, 0.08, 0.54), new = TRUE, mar = c(0, 0, 0, 0))
  draw_report_table(table, "Calculated-power diagnostics")
}

write_sensitivity_report <- function(file, fit) {
  table <- report_sensitivity_table(fit)

  grDevices::pdf(file, width = 8.5, height = 11)
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::text(0.5, 0.96, "ANOVA sensitivity analysis report", font = 2, cex = 1.25)
  graphics::text(0.06, 0.91, paste("Between-subjects factors:", report_factor_text(fit$design, "between")), adj = 0, cex = 0.8)
  graphics::text(0.06, 0.875, paste("Within-subjects factors:", report_factor_text(fit$design, "within")), adj = 0, cex = 0.8)
  graphics::text(0.06, 0.84, paste("Tested term:", fit$term), adj = 0, cex = 0.8)
  graphics::text(0.06, 0.805, report_fixed_n_text(fit), adj = 0, cex = 0.8)
  graphics::text(
    0.06, 0.77,
    sprintf("Target power: %.3f    Alpha: %.3f", fit$power, fit$alpha),
    adj = 0, cex = 0.8
  )
  graphics::text(
    0.06, 0.735,
    sprintf(
      "Effect-size search: [%.6g, %.6g]    Tolerance: %.6g",
      fit$pes_min, fit$pes_max, fit$pes_tol
    ),
    adj = 0, cex = 0.8
  )
  graphics::text(0.06, 0.70, report_convention_text(fit), adj = 0, cex = 0.8)
  sensitivity_text <- if (is.finite(fit$pes_needed)) {
    sprintf(
      "Minimum detectable partial eta squared: %.6f%s",
      fit$pes_needed,
      if (isTRUE(fit$converged)) "" else " (search did not converge)"
    )
  } else {
    "Target power was not reached within the effect-size search interval."
  }
  graphics::text(0.06, 0.665, sensitivity_text, adj = 0, font = 2, cex = 0.8)
  draw_report_footer()
  graphics::par(
    fig = c(0.11, 0.91, 0.45, 0.64),
    new = TRUE,
    mar = c(2.7, 3.5, 0.5, 0.5),
    mgp = c(1.7, 0.5, 0)
  )
  draw_sensitivity_curve(fit, legend_cex = 0.8)
  graphics::par(fig = c(0.06, 0.94, 0.04, 0.38), new = TRUE, mar = c(0, 0, 0, 0))
  draw_report_table(table, "Effect-size search results")
}

base64_encode <- function(bytes) {
  if (!length(bytes)) {
    return("")
  }
  alphabet <- strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/", "", fixed = TRUE)[[1]]
  padding <- (3L - length(bytes) %% 3L) %% 3L
  triplets <- matrix(
    c(as.integer(bytes), rep.int(0L, padding)),
    ncol = 3L,
    byrow = TRUE
  )
  indices <- cbind(
    bitwShiftR(triplets[, 1L], 2L),
    bitwOr(bitwShiftL(bitwAnd(triplets[, 1L], 3L), 4L), bitwShiftR(triplets[, 2L], 4L)),
    bitwOr(bitwShiftL(bitwAnd(triplets[, 2L], 15L), 2L), bitwShiftR(triplets[, 3L], 6L)),
    bitwAnd(triplets[, 3L], 63L)
  ) + 1L
  encoded <- vapply(
    seq_len(nrow(indices)),
    function(i) paste0(alphabet[indices[i, ]], collapse = ""),
    character(1)
  )
  if (padding == 1L) {
    encoded[[length(encoded)]] <- paste0(substr(encoded[[length(encoded)]], 1L, 3L), "=")
  } else if (padding == 2L) {
    encoded[[length(encoded)]] <- paste0(substr(encoded[[length(encoded)]], 1L, 2L), "==")
  }
  paste0(encoded, collapse = "")
}

ui <- fluidPage(
  tags$head(
    tags$title("anovapowersim Shiny app"),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('downloadPdf', function(message) {
        const binary = atob(message.content);
        const bytes = new Uint8Array(binary.length);
        for (let i = 0; i < binary.length; i++) {
          bytes[i] = binary.charCodeAt(i);
        }
        const url = URL.createObjectURL(new Blob([bytes], { type: 'application/pdf' }));
        const link = document.createElement('a');
        link.href = url;
        link.download = message.filename;
        document.body.appendChild(link);
        link.click();
        link.remove();
        setTimeout(function() { URL.revokeObjectURL(url); }, 0);
      });
    ")),
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
      .analysis-tabs .nav-tabs {
        border-bottom: 1px solid var(--line);
        margin: 0 0 16px;
      }
      .analysis-tabs .nav-tabs > li > a {
        color: var(--muted);
        padding: 8px 10px;
      }
      .analysis-tabs .nav-tabs > li.active > a,
      .analysis-tabs .nav-tabs > li.active > a:focus,
      .analysis-tabs .nav-tabs > li.active > a:hover {
        color: var(--accent-dark);
        font-weight: 650;
      }
      .analysis-help {
        color: var(--muted);
        font-size: 13px;
        margin: -4px 0 12px;
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
      .report-download {
        margin: 0 0 12px;
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
        p(class = "subtitle", "Balanced factorial ANOVA power analyses")
      ),
      tags$a(
        href = "https://shaheedazaad.github.io/anovapowersim/articles/calculated-power.html",
        target = "_blank",
        rel = "noopener noreferrer",
        "Calculated power guide"
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
        " R package. Browser runs use analytic calculations for required sample size, achieved power and sensitivity. ",
        "The R call below each analysis shows the equivalent package function for reproducible workflows."
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
          h2("Analysis"),
          selectInput("term", "Term", choices = character()),
          div(
            class = "analysis-tabs",
            tabsetPanel(
              id = "analysis_type",
              type = "tabs",
              selected = "sample_size",
              tabPanel(
                "Sample size",
                value = "sample_size",
                p(
                  class = "analysis-help",
                  "Find the sample size required for a target effect size and power."
                ),
                div(
                  class = "inline-grid",
                  numericInput("target_pes", "Target partial eta squared", value = 0.03, min = 0.001, max = 0.999, step = 0.001),
                  numericInput("power", "Target power", value = 0.90, min = 0.001, max = 0.999, step = 0.01),
                  numericInput("alpha", "Alpha", value = 0.05, min = 0.001, max = 0.999, step = 0.001),
                  numericInput("epsilon", "Nonsphericity correction", value = power_n_defaults$epsilon, min = 0.001, max = 1, step = 0.01),
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
                  actionButton("run", "Calculate sample size", icon = icon("play"), class = "btn-primary"),
                  uiOutput("run_status")
                )
              ),
              tabPanel(
                "Achieved power",
                value = "achieved",
                p(
                  class = "analysis-help",
                  "Calculate power for a fixed sample size and effect size."
                ),
                div(
                  class = "inline-grid",
                  numericInput("achieved_target_pes", "Partial eta squared", value = 0.08, min = 0.001, max = 0.999, step = 0.001),
                  numericInput("achieved_n", "n per between-subject cell", value = 30, min = 2, step = 1),
                  numericInput("achieved_alpha", "Alpha", value = power_effect_defaults$alpha, min = 0.001, max = 0.999, step = 0.001),
                  numericInput("achieved_epsilon", "Nonsphericity correction", value = power_effect_defaults$epsilon, min = 0.001, max = 1, step = 0.01)
                ),
                p(
                  class = "analysis-help",
                  "For a within-only design, n is the total sample size."
                ),
                tags$details(
                  class = "advanced-settings",
                  tags$summary("Advanced settings"),
                  checkboxInput("achieved_gpower", "G*Power convention", value = FALSE)
                ),
                div(
                  class = "run-row",
                  actionButton("run_achieved", "Calculate achieved power", icon = icon("play"), class = "btn-primary"),
                  uiOutput("achieved_run_status")
                )
              ),
              tabPanel(
                "Sensitivity",
                value = "sensitivity",
                p(
                  class = "analysis-help",
                  "Find the minimum detectable partial eta squared at a fixed sample size."
                ),
                div(
                  class = "inline-grid",
                  numericInput("sensitivity_n", "n per between-subject cell", value = 30, min = 2, step = 1),
                  numericInput("sensitivity_power", "Target power", value = power_effect_defaults$power, min = 0.001, max = 0.999, step = 0.01),
                  numericInput("sensitivity_alpha", "Alpha", value = power_effect_defaults$alpha, min = 0.001, max = 0.999, step = 0.001),
                  numericInput("sensitivity_epsilon", "Nonsphericity correction", value = power_effect_defaults$epsilon, min = 0.001, max = 1, step = 0.01)
                ),
                p(
                  class = "analysis-help",
                  "For a within-only design, n is the total sample size."
                ),
                tags$details(
                  class = "advanced-settings",
                  tags$summary("Search settings"),
                  div(
                    class = "inline-grid",
                    numericInput("pes_min", "Minimum partial eta squared", value = power_effect_defaults$pes_min, min = 0.000001, max = 0.999, step = 0.001),
                    numericInput("pes_max", "Maximum partial eta squared", value = power_effect_defaults$pes_max, min = 0.001, max = 0.999999, step = 0.01),
                    numericInput("pes_tol", "Search tolerance", value = power_effect_defaults$pes_tol, min = 0.000001, step = 0.001)
                  ),
                  checkboxInput("sensitivity_gpower", "G*Power convention", value = FALSE)
                ),
                div(
                  class = "run-row",
                  actionButton("run_sensitivity", "Run sensitivity analysis", icon = icon("play"), class = "btn-primary"),
                  uiOutput("sensitivity_run_status")
                )
              )
            )
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
          uiOutput("report_download"),
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
  running <- reactiveVal(NULL)
  result <- reactiveVal(NULL)
  result_type <- reactiveVal(NULL)
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
    current <- preserve_factor_inputs(factors(), input)
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
    current <- preserve_factor_inputs(factors(), input)
    factors(current[current$id != id, , drop = FALSE])
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

  design_call_args <- reactive({
    problem <- design_problem()
    validate(need(is.null(problem), problem))
    df <- current_design()
    req(input$term)
    list(
      between = named_counts(df, "between"),
      within = named_counts(df, "within"),
      term = input$term
    )
  })

  power_args <- reactive({
    c(
      design_call_args(),
      list(
        target_pes = input$target_pes,
        power = input$power,
        alpha = input$alpha,
        n_start = if (isTRUE(input$use_n_start)) as.integer(input$n_start) else NULL,
        n_max = as.integer(input$n_max),
        gpower = isTRUE(input$gpower),
        epsilon = input$epsilon
      )
    )
  })

  achieved_args <- reactive({
    c(
      design_call_args(),
      list(
        target_pes = input$achieved_target_pes,
        n = as.integer(input$achieved_n),
        alpha = input$achieved_alpha,
        gpower = isTRUE(input$achieved_gpower),
        epsilon = input$achieved_epsilon
      )
    )
  })

  sensitivity_args <- reactive({
    c(
      design_call_args(),
      list(
        n = as.integer(input$sensitivity_n),
        power = input$sensitivity_power,
        alpha = input$sensitivity_alpha,
        pes_min = input$pes_min,
        pes_max = input$pes_max,
        pes_tol = input$pes_tol,
        gpower = isTRUE(input$sensitivity_gpower),
        epsilon = input$sensitivity_epsilon
      )
    )
  })

  output$call <- renderText({
    tryCatch(
      switch(
        input$analysis_type %||% "sample_size",
        sample_size = format_call(
          "power_n_calc",
          power_n_call_args(power_args())
        ),
        achieved = format_call(
          "power_achieved_calc",
          power_achieved_call_args(achieved_args())
        ),
        sensitivity = format_call(
          "power_sensitivity_calc",
          power_sensitivity_call_args(sensitivity_args())
        ),
        ""
      ),
      error = function(e) ""
    )
  })

  run_analysis <- function(kind, backend, args) {
    # R function arguments are lazy. Force reactive-derived `args` while this
    # helper is still running inside the button observer; otherwise the first
    # evaluation would occur in `onFlushed()`, outside a reactive context.
    force(kind)
    force(backend)
    force(args)

    captured_warnings <- character()
    run_warnings(character())
    result(NULL)
    result_type(kind)
    running(kind)
    set_progress_state(list(
      active = TRUE,
      lines = character()
    ))

    session$onFlushed(function() {
      on.exit({
        running(NULL)
        set_progress_state(list(
          active = FALSE,
          lines = character()
        ))
      }, add = TRUE)

      fit <- withCallingHandlers(
        tryCatch(
          do.call(backend, args),
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
  }

  observeEvent(input$run, {
    run_analysis(
      kind = "sample_size",
      backend = power_n_backend,
      args = power_n_calc_args(power_args())
    )
  })

  observeEvent(input$run_achieved, {
    run_analysis(
      kind = "achieved",
      backend = power_achieved_backend,
      args = achieved_args()
    )
  })

  observeEvent(input$run_sensitivity, {
    run_analysis(
      kind = "sensitivity",
      backend = power_sensitivity_backend,
      args = sensitivity_args()
    )
  })

  output$run_status <- renderUI({
    if (identical(running(), "sample_size")) {
      span(class = "text-muted", "Calculating...")
    } else {
      NULL
    }
  })

  output$achieved_run_status <- renderUI({
    if (identical(running(), "achieved")) {
      span(class = "text-muted", "Calculating...")
    } else {
      NULL
    }
  })

  output$sensitivity_run_status <- renderUI({
    if (identical(running(), "sensitivity")) {
      span(class = "text-muted", "Calculating...")
    } else {
      NULL
    }
  })

  observeEvent(input$analysis_type, {
    result(NULL)
    result_type(NULL)
    run_warnings(character())
  }, ignoreInit = TRUE)

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
      empty_copy <- switch(
        input$analysis_type %||% "sample_size",
        sample_size = list(
          "Ready to calculate sample size",
          "Enter the design, effect size and target power, then click Calculate sample size."
        ),
        achieved = list(
          "Ready to calculate achieved power",
          "Enter a fixed sample size and effect size, then click Calculate achieved power."
        ),
        sensitivity = list(
          "Ready to run a sensitivity analysis",
          "Enter a fixed sample size and target power, then click Run sensitivity analysis."
        )
      )
      return(div(
        class = "empty-results",
        h3(empty_copy[[1L]]),
        p(empty_copy[[2L]])
      ))
    }

    kind <- term_kind(fit$term)
    alpha_text <- format_decimal(fit$alpha)
    between_cells <- fit$design$n_between_cells %||% 1L

    if (identical(result_type(), "sample_size")) {
      power_text <- format_decimal(fit$power)
      pes_text <- format_decimal(fit$target_pes)
      n_per_cell_text <- format_count(fit$n_needed)
      total_n_text <- format_count(fit$total_n_needed)
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

      return(div(
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
      ))
    }

    if (identical(result_type(), "achieved")) {
      n_detail <- if (between_cells > 1L) {
        sprintf(
          "%s participants per between-subjects cell; %s participants in total.",
          format_count(fit$n),
          format_count(fit$total_n)
        )
      } else {
        sprintf("%s participants in total.", format_count(fit$total_n))
      }
      return(div(
        class = "result-statement",
        p(
          class = "result-sentence",
          "Achieved power is ",
          strong(format_decimal(fit$achieved_power)),
          " for the ",
          strong(fit$term),
          " ",
          kind,
          "."
        ),
        p(
          class = "result-detail",
          sprintf(
            "%s Alpha = %s and partial eta squared = %s.",
            n_detail,
            alpha_text,
            format_decimal(fit$target_pes)
          )
        )
      ))
    }

    pes_reached <- is.finite(fit$pes_needed)
    div(
      class = if (pes_reached) "result-statement" else "result-statement not-reached",
      if (pes_reached) {
        p(
          class = "result-sentence",
          "The minimum detectable partial eta squared is ",
          strong(format_decimal(fit$pes_needed, digits = 4)),
          " for the ",
          strong(fit$term),
          " ",
          kind,
          "."
        )
      } else {
        p(
          class = "result-sentence",
          "The target power of ",
          strong(format_decimal(fit$power)),
          " was not reached within the effect-size search interval."
        )
      },
      p(
        class = "result-detail",
        sprintf(
          "Based on %s total participants, target power = %s, alpha = %s, and a final bracket width of %s.",
          format_count(fit$total_n),
          format_decimal(fit$power),
          alpha_text,
          if (is.finite(fit$bracket_width)) format(fit$bracket_width, digits = 4) else "not available"
        )
      )
    )
  })

  output$result_view_ui <- renderUI({
    req(result())
    if (identical(result_type(), "achieved")) {
      return(div(
        class = "result-view table-scroll",
        tableOutput("results_table")
      ))
    }
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

  output$report_download <- renderUI({
    req(result())
    label <- switch(
      result_type(),
      sample_size = "Download sample-size PDF report",
      achieved = "Download achieved-power PDF report",
      sensitivity = "Download sensitivity PDF report",
      "Download PDF report"
    )
    div(
      class = "report-download",
      actionButton("download_report", label, icon = icon("file-pdf"))
    )
  })

  observeEvent(input$download_report, {
    fit <- result()
    req(fit)
    analysis <- result_type()
    req(analysis %in% c("sample_size", "achieved", "sensitivity"))
    report_file <- tempfile(fileext = ".pdf")
    on.exit(unlink(report_file), add = TRUE)
    tryCatch({
      report_writer <- switch(
        analysis,
        sample_size = write_power_report,
        achieved = write_achieved_report,
        sensitivity = write_sensitivity_report
      )
      report_writer(report_file, fit)
      report_slug <- switch(
        analysis,
        sample_size = "sample-size",
        achieved = "achieved-power",
        sensitivity = "sensitivity"
      )
      session$sendCustomMessage("downloadPdf", list(
        filename = sprintf(
          "anovapowersim-%s-report-%s.pdf",
          report_slug,
          format(Sys.Date(), "%Y-%m-%d")
        ),
        content = base64_encode(readBin(report_file, what = "raw", n = file.info(report_file)$size))
      ))
    }, error = function(e) {
      showNotification(
        paste("Could not create the PDF report:", conditionMessage(e)),
        type = "error",
        duration = 10
      )
    })
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
    if (identical(result_type(), "sensitivity")) {
      draw_sensitivity_curve(fit)
    } else {
      draw_power_curve(fit)
    }
  })

  output$results_table <- renderTable({
    fit <- result()
    req(fit)
    results <- as.data.frame(fit$results)
    simulation_fields <- c("n_sims", "valid_sims", "failed_sims", "power_sim")
    results <- results[, setdiff(names(results), simulation_fields), drop = FALSE]
    if ("target_pes" %in% names(results)) {
      results$target_pes <- sprintf("%.6f", results$target_pes)
    }
    results
  }, striped = TRUE, bordered = TRUE, digits = 3)
}

shinyApp(ui, server)
