shinyServer(function(input, output, session) {

  dataStore    <- reactiveVal(list())
  dataMetaInfo <- reactiveVal(list())
  genCounter   <- reactiveVal(0)
  previewState <- reactiveValues()
  fitStateMap  <- reactiveValues()

  makeSafeId <- function(dname) {
    paste0("ds_", digest::digest(dname, algo = "md5"))
  }

  output$gen_params_ui <- renderUI({
    req(input$gen_distr)
    distr    <- getDistrInfo(input$gen_distr)
    defaults <- distr$param.defaults

    lapply(distr$param.names, function(nm) {
      numericInput(paste0('genparam_', input$gen_distr, '_', nm), 
                   label = nm,
                   value = defaults[[nm]], step = 0.1)
    })
  })

  observeEvent(input$gen_btn, {
    req(input$gen_distr, input$gen_n)
    distr  <- getDistrInfo(input$gen_distr)
    pnames <- distr$param.names

    params <- lapply(pnames, function(nm) input[[paste0('genparam_', input$gen_distr, '_', nm)]])
    names(params) <- pnames

    if (any(sapply(params, is.null))) {
      showNotification('Distribution parameters missing -- try again.', type = 'error')
      return()
    }

    if (!is.na(as.integer(input$gen_seed))) {
      set_seed(as.integer(input$gen_seed)) 
    } else {
      set_seed(NULL)
    }

    xfunc <- sprintf('x%s', distr$name)
    args  <- c(list(n = input$gen_n), params)
          
    x <- tryCatch(do.call(xfunc, args), error = function(e) NULL)

    if (is.null(x)) {
      showNotification("Error generating data series.", type = "error")
      return()
    }

    count <- genCounter() + 1
    genCounter(count)
    dataset_name <- sprintf("Synthetic Data Set %d: %s (n=%d)", 
                             count, input$gen_distr, input$gen_n)
    dataset_meta <- list(idx = count,
                         distr = input$gen_distr,
                         n = input$gen_n,
                         params = params)

    current_list <- dataStore()
    current_list[[dataset_name]] <- x
    dataStore(current_list)
    
    metainfo <- dataMetaInfo()
    metainfo[[count]] <- dataset_meta
    dataMetaInfo(metainfo)
    
    fitStateMap[[dataset_name]]  <- FALSE
    previewState[[dataset_name]] <- FALSE
  })

  observeEvent(input$datafile, {
    req(input$datafile)
    filepath <- input$datafile$datapath

    df <- tryCatch({
      read.csv(filepath, header = input$header, sep = "", check.names = FALSE)
    }, error = function(e) {
      tryCatch(read.csv(filepath, header = input$header, check.names = FALSE), error = function(ev) NULL)
    })

    if (is.null(df) || ncol(df) == 0) {
      showNotification("Failed to parse file. Make sure file format is TSV, CSV, or Tab-delimited text.", type = "error")
      return()
    }

    current_list <- dataStore()
    num_cols <- 0
    for (col_name in names(df)) {
      vec <- suppressWarnings(as.numeric(df[[col_name]]))
      vec <- vec[!is.na(vec)]
      if (length(vec) >= MIN_OBS) {
        num_cols <- num_cols + 1
        entry_key <- sprintf("File: %s [%s]", input$datafile$name, col_name)
        current_list[[entry_key]] <- vec
        fitStateMap[[entry_key]]  <- FALSE
        previewState[[entry_key]] <- FALSE
      }
    }
    
    if (num_cols == 0) {
      showNotification("No numeric columns with adequate observations were found.", type = "warning")
    } else {
      dataStore(current_list)
    }
  })

  observeEvent(input$gen_reset, {
    dataStore(list())
    genCounter(0)
    showNotification("All generated data cleared.", type = "message")
  })

  output$copy_button_ui <- renderUI({
    ds <- dataStore()
    if (length(ds) == 0) return(NULL)
    tsv_text <- list2tsv(ds)
    rclipButton(
      inputId = "copy_btn",
      label = "Copy All Data",
      clipText = tsv_text,
      icon = icon("copy")
    )
  })

  output$dl_all_btn <- downloadHandler(
    filename = function() {
      paste0("all_generated_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv")
    },
    content = function(file) {
      ds <- dataStore()
      tsv_text <- formatDatasetsToTSV(ds)
      writeLines(tsv_text, file)
    }
  )

  getTopFits <- function(x, pcutoff, alpha, selected_distrs) {
    distrs <- if (is.null(selected_distrs) || "All applicable" %in% selected_distrs) NULL else selected_distrs
    fits <- tryCatch(
      distrEstimate(x = x, distribution = distrs, alpha = alpha, return.best = FALSE),
      error = function(e) NULL
    )
    if (is.null(fits) || length(fits) == 0) return(NULL)

    pvals      <- sapply(fits, function(r) r$p.value)
    res.sorted <- fits[order(pvals, decreasing = TRUE)]
    good       <- res.sorted[sapply(res.sorted, function(r) r$p.value >= pcutoff)]

    if (length(good) > 0) {
      good[seq_len(min(3, length(good)))]
    } else {
      res.sorted[seq_len(min(3, length(res.sorted)))]
    }
  }

  # --- HTML REPORT GENERATOR WITH EMBEDDED PLOTS ---
  output$dl_report_html <- downloadHandler(
    filename = function() {
      paste0("fitting_report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
    },
    content = function(file) {
      ds <- dataStore()
      if (length(ds) == 0) {
        writeLines("<html><body><h3>No datasets available to generate report.</h3></body></html>", file)
        return()
      }

      html_content <- c(
        "<!DOCTYPE html>",
        "<html><head><title>Fitting Results Report</title>",
        "<style>",
        "body { font-family: Arial, sans-serif; margin: 30px; background-color: #f9f9f9; }",
        ".card { background: white; padding: 25px; margin-bottom: 25px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }",
        "table { border-collapse: collapse; width: 100%; margin-top: 15px; }",
        "th, td { border: 1px solid #ddd; padding: 10px; text-align: left; }",
        "th { background-color: #007bff; color: white; }",
        ".plot-img { max-width: 100%; height: auto; display: block; margin: 15px 0; border: 1px solid #e0e0e0; border-radius: 4px; }",
        "</style></head><body>",
        "<h1>Distribution Fitting Report</h1>",
        sprintf("<p><em>Generated on: %s</em></p><hr>", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
      )

      for (dname in names(ds)) {
        x <- ds[[dname]]
        sid <- makeSafeId(dname)
        is_fitted <- isTRUE(fitStateMap[[dname]])
        
        pcutoff <- input[[paste0("pcutoff_", sid)]] %||% 0.1
        alpha   <- input[[paste0("alpha_", sid)]] %||% 0.05
        fit_distrs <- input[[paste0("fit_distrs_", sid)]] %||% "All applicable"

        top_fits <- if (is_fitted) getTopFits(x, pcutoff, alpha, fit_distrs) else NULL

        # Render Plot to base64 string
        tmp_img <- tempfile(fileext = ".png")
        png(tmp_img, width = 800, height = 500, res = 100)
        plotDistr(x, top.results = top_fits, is.fitted = is_fitted, title = dname)
        dev.off()

        img_b64 <- base64encode(tmp_img)
        unlink(tmp_img)

        # Build Card
        html_content <- c(html_content, sprintf("<div class='card'><h2>Dataset: %s</h2>", dname))
        html_content <- c(html_content, sprintf("<p><strong>Observations:</strong> %d | <strong>Mean:</strong> %.4f | <strong>SD:</strong> %.4f</p>", length(x), mean(x), sd(x)))
        
        # Embedded Plot Image
        html_content <- c(html_content, sprintf("<img class='plot-img' src='data:image/png;base64,%s' />", img_b64))

        if (is_fitted) {
          if (!is.null(top_fits) && length(top_fits) > 0) {
            html_content <- c(html_content, "<h3>Top Distribution Fits</h3><table><tr><th>Rank</th><th>Distribution</th><th>P-Value</th><th>Chi-Square</th><th>Parameters</th></tr>")
            for (k in seq_along(top_fits)) {
              r <- top_fits[[k]]
              params <- as.list(r$parameters)
              param_str <- paste(sapply(names(params), function(p) sprintf("%s=%.4f", p, params[[p]])), collapse = ", ")
              pval_str <- if (r$p.value < 0.0001) sprintf("%.4e", r$p.value) else sprintf("%.4f", r$p.value)
              html_content <- c(html_content, sprintf("<tr><td>#%d</td><td>%s</td><td>%s</td><td>%.2f</td><td>%s</td></tr>",
                                                     k, rev(r$distr.name)[1], pval_str, r$chisq.statistics, param_str))
            }
            html_content <- c(html_content, "</table>")
          } else {
            html_content <- c(html_content, "<p><em>No valid distributions passed the fitting criteria.</em></p>")
          }
        } else {
          html_content <- c(html_content, "<p><em>No fitting detail available.</em></p>")
        }
        html_content <- c(html_content, "</div>")
      }

      html_content <- c(html_content, "</body></html>")
      writeLines(html_content, file)
    }
  )

  # --- UI RENDERER & DYNAMIC DATA BINDINGS ---
  output$dynamicDataRows <- renderUI({
    ds <- dataStore()
    dsMeta <- dataMetaInfo()
    if (length(ds) == 0) {
      return(wellPanel(helpText("No datasets generated or uploaded yet. Click 'Generate Numbers' or upload a file.")))
    }

    items <- lapply(names(ds), function(dname) {
      safe_id <- makeSafeId(dname)
      show_preview <- isTRUE(previewState[[dname]])

      # Register outputs dynamically per dataset
      local({
        target_name <- dname
        sid         <- safe_id

        output[[paste0("combined_plot_", sid)]] <- renderPlot({
          current_ds <- dataStore()
          vec <- current_ds[[target_name]]
          req(!is.null(vec))
          is_fitted <- isTRUE(fitStateMap[[target_name]])
          
          fit_distrs <- input[[paste0("fit_distrs_", sid)]]
          alpha      <- input[[paste0("alpha_", sid)]] %||% 0.05
          pcutoff    <- input[[paste0("pcutoff_", sid)]] %||% 0.1

          top_fits <- if (is_fitted) getTopFits(vec, pcutoff, alpha, fit_distrs) else NULL
          plotDistr(vec, top.results = top_fits, is.fitted = is_fitted, title = target_name)
        })

        output[[paste0("details_", sid)]] <- renderText({
          is_fitted <- isTRUE(fitStateMap[[target_name]])
          if (!is_fitted) return("no fitting detail available")
          
          current_ds <- dataStore()
          vec <- current_ds[[target_name]]
          if (is.null(vec)) return("")

          fit_distrs <- input[[paste0("fit_distrs_", sid)]]
          alpha      <- input[[paste0("alpha_", sid)]] %||% 0.05
          pcutoff    <- input[[paste0("pcutoff_", sid)]] %||% 0.1

          top_fits <- getTopFits(vec, pcutoff, alpha, fit_distrs)
          if (is.null(top_fits) || length(top_fits) == 0) {
            return("No valid distributions fit for this data.")
          }
          
          out <- sapply(seq_along(top_fits), function(k) {
            r <- top_fits[[k]]
            pval <- r$p.value
            pval_str <- if (pval < 0.0001) sprintf("%.4e", pval) else sprintf("%.4f", pval)
            params <- as.list(r$parameters)
            param_lines <- sapply(names(params), function(pname) {
              sprintf("  %-10s = %.4f", pname, params[[pname]])
            })
            param_str <- paste(param_lines, collapse = "\n")
            sprintf("#%d Distribution: %s\n   P-Value: %s\n   ChiSq:   %.2f\nParameters:\n%s",
                    k, rev(r$distr.name)[1], pval_str, r$chisq.statistics, param_str)
          })
          paste(out, collapse = "\n-----------------------------\n")
        })

      })

      ## render data row
      wellPanel(
        fluidRow(
          column(8, h4(dname)),
          column(4, class = "text-right",
                 align = 'left',
                 actionButton(paste0("btn_rem_", safe_id), 
                              "Remove Data", 
                              class = "btn-sm btn-danger", icon = icon("trash"))),
        ),
        hr(),
        
        fluidRow(
          column(6, 
                 plotOutput(paste0("combined_plot_", safe_id), height = "380px")),
          
          column(3, 
                 h5("Fit Details (Top 3)"),
                 verbatimTextOutput(paste0("details_", safe_id))),
          column(3,
                 fluidRow(
                          selectInput(paste0("fit_distrs_", safe_id),
                                      "Distributions to test",
                                      choices = c("All applicable", DISTR_CHOICES), 
                                      selected = "All applicable", 
                                      multiple = TRUE),
                          numericInput(paste0("alpha_", safe_id), 
                                       "Chi-sq test alpha", 
                                       value = 0.05, 
                                       min = 0.001, 
                                       max = 0.5, 
                                       step = 0.01),
                          numericInput(paste0("pcutoff_", safe_id), 
                                       "P-value cutoff",
                                       value = 0.1, 
                                       min = 0, 
                                       max = 0.999, 
                                       step = 0.01),
                          actionButton(paste0("btn_fit_", safe_id), 
                                       "Fit Data", 
                                       class = "btn-success btn-block", 
                                       icon = icon("chart-line"))
                          )
                )

         ) ## end of data row
      )    ## end of wellPanel
    })
    
    do.call(tagList, items)
  })

  # --- GLOBAL EVENT DELEGATOR FOR PREVIEW, FIT & REMOVE BUTTONS ---
  observe({
    ds <- dataStore()
    if (length(ds) == 0) return()

    for (dname in names(ds)) {
      local({
        target_name <- dname
        sid         <- makeSafeId(target_name)
        
        btn_rem_id  <- paste0("btn_rem_", sid)
        btn_fit_id  <- paste0("btn_fit_", sid)

        observeEvent(input[[btn_fit_id]], {
          fitStateMap[[target_name]] <- TRUE
        }, ignoreInit = TRUE, autoDestroy = TRUE)

        observeEvent(input[[btn_rem_id]], {
          current_ds <- dataStore()
          current_ds[[target_name]] <- NULL
          dataStore(current_ds)
          showNotification(sprintf("Removed dataset: %s", target_name), type = "message")
        }, ignoreInit = TRUE, autoDestroy = TRUE)
      })
    }
  })

})