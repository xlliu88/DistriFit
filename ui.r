shinyUI(fluidPage(

  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('copyToClipboard', function(text) {
        if (navigator.clipboard && window.isSecureContext) {
          navigator.clipboard.writeText(text);
        } else {
          var textArea = document.createElement('textarea');
          textArea.value = text;
          textArea.style.position = 'fixed';
          textArea.style.left = '-999999px';
          document.body.appendChild(textArea);
          textArea.focus();
          textArea.select();
          try {
            document.execCommand('copy');
          } catch (err) {
            console.error('Fallback copy failed: ', err);
          }
          document.body.removeChild(textArea);
        }
      });
    "))
  ),
  rclipboardSetup(),
  titlePanel("Pseudo Random Number Generator | Distribution Fitting"),

  sidebarLayout(

    sidebarPanel(
      width = 3,

      tabsetPanel(
        id = "input_tabs",

        tabPanel("Generate Data",
          br(),
          selectInput("gen_distr", 
                      "Distribution to sample from",
                      choices = DISTR_CHOICES, 
                      selected = "unif"),
          uiOutput("gen_params_ui"),
          numericInput("gen_n", 
                       "Number of observations", 
                       value = 100, 
                       min = 30, 
                       step = 1),
          textInput("gen_seed", 
                    "Random seed (optional)", 
                    value = 'Auto'),
          br(),
          actionButton("gen_btn", 
                       "Generate Numbers", 
                       class = "btn-primary", 
                       icon = icon("dice")),
          actionButton("gen_reset", 
                       "Clear All", 
                       class = 'btn-danger', 
                       icon = icon("rotate-left"))
        ),
        
        tabPanel("Upload File",
          br(),
          fileInput("datafile", 
                    "Upload CSV / TSV / TXT file", 
                    accept = c(".csv", ".tsv", ".txt")),
          checkboxInput("header", 
                        "File has header row", 
                        value = FALSE)
        )

      ),

      hr(),
      h4("Fitting options"),

      selectInput("fit_distrs",
                  "Distributions to test",
                  choices = c("All applicable", DISTR_CHOICES), 
                  selected = "All applicable", 
                  multiple = TRUE),

      numericInput("alpha", 
                   "Chi-sq test alpha", 
                   value = 0.05, 
                   min = 0.001, 
                   max = 0.5, 
                   step = 0.01),

      numericInput("pcutoff", 
                   "P-value cutoff",
                   value = 0.1, 
                   min = 0, 
                   max = 0.999, 
                   step = 0.01),

      br(),
      actionButton("fit_btn", 
                   "Goodness of fit", 
                   class = "btn-success", 
                   icon = icon("chart-line")),
      hr(),
      
      h5("Global Data Actions"),
      uiOutput("copy_button_ui"),
      downloadButton("dl_all_btn", 
                     "Save All Data", 
                     class = "btn-block btn-default"),
      
      downloadButton("dl_report_html", 
                     "Save Fitting Report (HTML)", 
                     class = "btn-block btn-info", 
                     icon = icon("file-code"))
      
    ),

    mainPanel(
      width = 9,
      uiOutput("dynamicDataRows")
    )
  )
))