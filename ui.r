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
                      selected = "norm"),
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