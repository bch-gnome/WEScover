# load packages
library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(shinyjs)
library(shinyBS)
library(reshape2)
library(RColorBrewer)
library(fst)
library(data.table)
library(wiggleplotr)
library(patchwork)
library(ggpubr)
library(dplyr)
library(corrplot)
setDTthreads(18)

# laod functions
source("function.R")

# load data
gene_symbol <- read.fst("data/gene_symbol.fst")
gnomad_exome <- read.fst("data/gnomad_exome.r2.1.fst")

gtrM <- read.fst("data/gtrM.fst")
gpt_tP_tG <- read.fst("data/gpt_tP_tG.fst")
ccds2ens <- readRDS("data/ccds_ens_map.rds")
# using comprehensive version of ORPHA data that includes genes not mapped to gpt
ORPHA <- read.fst("data/new_ORPHA.fst")
per_gene_summary<-read.fst("data/summary.long.fst")
kgp.map<-read.table("data/integrated_call_samples_v3.20130502.ALL.panel", header = T, sep = "\t", as.is = T, col.names = c("ID","sub_pop","pop","gender"))

load("data/tx_data.RData")

# format data
gpt_tP_tG$GTR_accession <- as.character(gpt_tP_tG$GTR_accession)
gpt_tP_tG$test_name <- as.character(gpt_tP_tG$test_name)
gpt_tP_tG$phenotype_name <- as.character(gpt_tP_tG$phenotype_name)
ORPHA$HPO_term_name <- as.character(ORPHA$HPO_term_name)
ORPHA$gene_symbol <- as.character(ORPHA$gene_symbol)

setDTthreads(18)

# Define UI ----
ui <- fluidPage(
  useShinyjs(),
  theme = shinytheme("flatly"),
  tags$head(tags$style(".modal-dialog{min-width:1200px}")),
  tags$head(includeHTML(("google_analytics.html"))),
  ## javascript required to clean the 'check' of detail buttons in main table
  tags$script("
      Shiny.addCustomMessageHandler('resetInputValue', function(variableName){
        Shiny.onInputChange(variableName, null);
      });
    "),
  tags$style(type="text/css", "body {padding-top: 80px;} .selectize-input {height: 45px;} .action-button {height:45px; width:100%;} .center {display: block; margin-left: auto; margin-right: auto}"),
  #js function to reset a button, variableName is the button name whose value we want to reset
  tags$script("Shiny.addCustomMessageHandler('reset_detail_button', function(detail_button){
                Shiny.onInputChange(detail_button, null);
                });
                "),
  navbarPage("WEScover", id="mainNav", windowTitle = "WEScover", position = "fixed-top", fluid = TRUE,
    tabPanel("Home",
     absolutePanel( width = "70%", left = "15%", right = "15%",
       wellPanel(
         em(h1("WEScover")),
         hr(),
         p(em('WEScover'), 'helps users to check whether genes of interest could be sufficiently covered in terms of breadth and depth by whole exome sequencing (WES). For each transcript, breadth of coverage data was calculated at various read depth from the ', 
           a("1000 Genomes Project (1KGP)", href = "http://www.internationalgenome.org/", target="_blank"), 
           '(N = 2,504). A user will be able to minimize the chance of false negatives by selecting a targeted gene panel test for the genes that WES cannot cover well.'),
         p('Breadth and depth of coverage for ', a(em('NOTCH1'), href = "http://gnomad.broadinstitute.org/gene/ENSG00000148400", target="_blank"),
           ' are illustrated below. For some of the exons, breadth of coverage seems to be sub-optimal that could result in false negative results with WES.'),
         tags$img(src=paste0("gnomAD_notch1.png?v=", as.numeric(Sys.time())), alt = "Coverage from gnomAD project for NOTCH1", style="width:650px;height:300px", class="center"),
         br(),
         p(em('WEScover'), ' provides detailed coverage information including difference in breadth of coverage between continent-level populatios.'),
         br(),
         tags$img(src=paste0("violin_notch1.png?v=", as.numeric(Sys.time())), alt = "Contintental population breath of coverage violin plot for CCDS43905.1/NOTCH1", style="width:650px;height:300px", class="center"),
         br(),
         p('Phenotype, genetic test names, or gene symbols can be used to retrieve coverage information in the query window. The output summary helps users to choose WES vs. targeted gene panel testing.')
       )
     )
    ),
    tabPanel("Query",
      sidebarLayout(position = "left",
      sidebarPanel(
              tags$h2("User input"),
              fluidRow(
                column(12, 
                       radioButtons("select_phen", "Select source of phenotype terms",
                                   c("Genetic Testing Registry (GTR)" = "GTR",
                                     "Human Phenotype Ontology (HPO)" = "HPO"),
                                   inline = FALSE)
                )
              ),
              fluidRow(
                column(8,
                       selectizeInput("phen",
                                      label="GTR Phenotype",
                                      choices = NULL,
                                      multiple = TRUE)
                ),
                column(4,
                       HTML("<label>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</label>"),
                       tags$br(),
                       actionButton("fGPT", "Filter")#, icon = icon("filter", lib = "glyphicon"))
                )
              ),
              # test input for HPO
              fluidRow(
                column(8,
                       selectizeInput("HPO",
                                      label="HPO Phenotype",
                                      choices = NULL,
                                      multiple = TRUE)
                ),
                column(4,
                       HTML("<label>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</label>"),
                       tags$br(),
                       actionButton("fGPT2", "Filter")#, icon = icon("filter", lib = "glyphicon"))
                       # actionButton("fGenes2", "Filter")#, icon = icon("filter", lib = "glyphicon"))
              )
              ),
              fluidRow(
                column(8,
                  selectizeInput("gpt",
                    label="GPT name",
                    choices = NULL,
                    multiple = TRUE)
                ),
                column(4,
                  HTML("<label>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</label>"),
                  tags$br(),
                  actionButton("fGenes", "Filter")#, icon = icon("filter", lib = "glyphicon"))
                )
              ),
              bsTooltip("gpt", title = "Select gene panel test(s) to find related genes", placement = "top", trigger = "hover"),
              bsTooltip("fGenes", title = "Works only after selecting phenotype(s)", placement = "right", trigger = "hover"),
              fluidRow(
                column(12, 
                  selectizeInput("gene_symbol", label = "Gene symbol",
                    choices = NULL,
                    multiple = TRUE,
                    options = list(
                      splitOn = I("(function() { return /[,; ]/; })()"),
                      create = I("function(input, callback){
                                         return {
                                         value: input,
                                         text: input
                                         };
                                         }")
                    )
                  )
                )
              ),
              fluidRow(
                column(12, 
                 selectInput("depth_of_coverage",
                             label = "Depth of coverage",
                             choices = c("5x", "10x", "15x", "20x", "25x", "30x", "50x", "100x"),
                             selected = "20x")
                             #width = '70%'),
                )
              ),
              fluidRow(
                column(12, 
                       radioButtons("assembly", "Human reference genome assembly version",
                                   c("GRCh37 (b37/hg19)" = "b37",
                                     "GRCh38 (hg38)" = "hg38"),
                                   selected = "b37",
                                   inline = TRUE)
                )
              ),
              fluidRow(
                column(2),
                column(5, actionButton("clear", "Clear inputs", class = "btn-secondary")),
                column(5, actionButton("update", "Submit query", class = "btn-primary"))
              )
              ),
        mainPanel(
          dataTableOutput('tableMain')
        )
      )
    ),
    tabPanel("Data",
      includeHTML("data.html"),
      hidden(
        numericInput( inputId = 'refresh_helper', label = 'refresh_helper', value = 0 ) 
      )
    ),
    tabPanel(title=HTML("</a></li><li><a href='https://bch-gnome.github.io/wescover_doc/master/' target='_blank'>Help")) ## workaround for linking help page
  )
)



# Define server logic ----
server <- function(input, output, session) {
  
  observeEvent(input$select_phen, {
    if (input$select_phen == "GTR") {
      updateSelectizeInput(session, 'phen', choices = sort(unique(gpt_tP_tG$phenotype_name)), server = TRUE)
      shinyjs::enable("phen")
      shinyjs::disable("HPO")
    }
    if (input$select_phen == "HPO") {
      updateSelectizeInput(session, 'HPO', choices = sort(unique(ORPHA$HPO_term_name)), server = TRUE)
      shinyjs::enable("HPO")
      shinyjs::disable("phen")
    }
  })
  
  # fill with default values
  updateSelectizeInput(session, 'HPO', choices = sort(unique(ORPHA$HPO_term_name)), server = TRUE)
  updateSelectizeInput(session, 'phen', choices = sort(unique(gpt_tP_tG$phenotype_name)), server = TRUE)
  updateSelectizeInput(session, 'gene_symbol', choices = sort(unique(gpt_tP_tG$gene_symbol)), server = TRUE)
  # updateSelectizeInput(session, 'gene_symbol', choices = sort(unique(genes_by_ccds_id$gene_symbol)), server = TRUE)
  updateSelectizeInput(session, 'gpt', choices = sort(unique(gpt_tP_tG$test_name)), server = TRUE)
  
  # if a gene symbol is provided by url
  observe({
    query <- parseQueryString(session$clientData$url_search)
    message("query -> ", query[['gene']])
    if (!is.null(query[['gene']]) & length(input$gene_symbol) == 0) {
      geneS <- strsplit(query[['gene']], ",")[[1]]
      updateNavbarPage(session, "mainNav", "Query")

      updateSelectizeInput(session, 'gene_symbol', choices = sort(unique(gpt_tP_tG$gene_symbol)), selected = query[['gene']], server = TRUE)
      # updateSelectizeInput(session, 'gene_symbol', choices = sort(unique(genes_by_ccds_id$gene_symbol)), selected = query[['gene']], server = TRUE)

    }
    if(length(input$gene_symbol) != 0 & input$refresh_helper == 0) {


    #}
    #if(length(input$gene_symbol) != 0 & input$refresh_helper == 0) {
      message("update! ", input$refresh_helper, " ", myValue$inc)

      updateNumericInput( session = session, inputId = 'refresh_helper', value = input$refresh_helper + 1 )
      # updateSelectizeInput(session, 'gene_symbol', choices = sort(unique(genes_by_ccds_id$gene_symbol)), selected = input$gene_symbol, server = TRUE)
      updateSelectizeInput(session, 'gene_symbol', choices = sort(unique(gpt_tP_tG$gene_symbol)), selected = input$gene_symbol, server = TRUE)
      if(input$refresh_helper %in% c(0,1)) {
        message("update from (RH) ", input$refresh_helper)
        updateNumericInput( session = session, inputId = 'refresh_helper', value = input$refresh_helper + 1 )
      }

    }
  })
  
  # if clear button is pushed
  observeEvent (input$clear,{
    updateSelectizeInput(session, 'HPO', choices = sort(unique(ORPHA$HPO_term_name)), server = TRUE)
    updateSelectizeInput(session, 'phen', choices = sort(unique(gpt_tP_tG$phenotype_name)), server = TRUE)

    # updateSelectizeInput(session, 'gene_symbol', 
    #                      choices = sort(unique(gpt_tP_tG$gene_symbol)), server = TRUE,
    #                      label = "Gene symbol")
    # updateSelectizeInput(session, 'gene_symbol', 
    #                      choices = gene_symbol$gene_symbol, server = TRUE,
    #                      label = "Gene symbol")
    updateSelectizeInput(session, 'gpt', choices = sort(unique(gpt_tP_tG$test_name)), server = TRUE,
                         label = "GPT name")
    # updateSelectizeInput(session, 'gene_symbol', choices = sort(unique(genes_by_ccds_id$gene_symbol)), server = TRUE, label = "Gene symbol")
    updateSelectizeInput(session, 'gene_symbol', choices = sort(unique(gpt_tP_tG$gene_symbol)), server = TRUE, label = "Gene symbol")

    # updateSelectizeInput(session, 'gpt', choices = sort(unique(gpt_tP_tG$test_name)), server = TRUE, label = "GPT name")
    updateSelectInput(session, "depth_of_coverage", choices = c("5x", "10x", "15x", "20x", "25x", "30x", "50x", "100x"), selected = "20x")
  })
  
  # if filter GPT button is pushed
  observeEvent (input$fGPT,{
    if (length(input$phen) != 0 & length(input$HPO) == 0) {
      listPhe <- unique(
        gpt_tP_tG$test_name[ gpt_tP_tG$GTR_accession %in% gpt_tP_tG$GTR_accession[ gpt_tP_tG$phenotype_name %in% as.character(input$phen) ] ])
      updateSelectizeInput(session, 'gpt', choices = sort(listPhe), server = TRUE, label = "GPT name (filtered)")
  
      listGenes <- unique(
        gpt_tP_tG$gene_symbol[ gpt_tP_tG$GTR_accession %in% gpt_tP_tG$GTR_accession[ gpt_tP_tG$phenotype_name %in% as.character(input$phen) ]])
      updateSelectizeInput(session, 'gene_symbol', choices = sort(listGenes), server = TRUE, label = "Gene symbol (filtered)")
    }
  })
  # if filter GPT2 button is pushed
  
  observeEvent (input$fGPT2,{
    if (length(input$phen) == 0 & length(input$HPO) != 0) {
      listPhe <- unique(
        ORPHA$test_name[ ORPHA$OrphaNumber_phen %in% ORPHA$OrphaNumber_phen[ ORPHA$HPO_term_name %in% as.character(input$HPO) ] ])
      updateSelectizeInput(session, 'gpt', choices = sort(listPhe), server = TRUE, label = "GPT name (filtered)")
  
      listGenes <- unique(
        ORPHA$gene_symbol[ ORPHA$OrphaNumber_phen %in% ORPHA$OrphaNumber_phen[ ORPHA$HPO_term_name %in% as.character(input$HPO) ]])
      updateSelectizeInput(session, 'gene_symbol', choices = sort(listGenes), server = TRUE, label = "Gene symbol (filtered)")
    }
  })

  
  # if filter genes button is pushed
  observeEvent (input$fGenes,{
    if (length(input$phen) != 0 & length(input$HPO) == 0) {
      listGenes <- gpt_tP_tG$gene_symbol[ gpt_tP_tG$test_name %in% as.character(input$gpt) ]
      updateSelectizeInput(session, 'gene_symbol', choices = sort(listGenes), server = TRUE, label = "Gene symbol (filtered)")
    }
    if (length(input$phen) == 0 & length(input$HPO) != 0) {
      listGenes <- ORPHA$gene_symbol[ ORPHA$test_name %in% as.character(input$gpt) ]
      updateSelectizeInput(session, 'gene_symbol', choices = sort(listGenes), server = TRUE, label = "Gene symbol (filtered)")
    }
  })

  # filter for genes in HPO
  # observeEvent (input$fGenes2,{
  #   listGenes <- ORPHA$gene_symbol[ ORPHA$HPO_term_name %in% as.character(input$HPO) ]
  #   updateSelectizeInput(session, 'gene_symbol', choices = sort(listGenes), server = TRUE, label = "Gene symbol (filtered)")
  # })
  
  # generator of buttons for the main table
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
  # list of global values
  myValue <- reactiveValues(summary_table = NA, GPT_table = NA,
    violon_population = NA, gene = "", ccds="", assembly="")
  
  # observer to capture detail buttons in the main table
  observeEvent(input$detail_button, {
      selectedRow <- as.numeric(strsplit(input$detail_button, "_")[[1]][2])
      geneS <- as.character(main_table()[selectedRow, 1])
      ccds <-  as.character(main_table()[selectedRow, 2])
      withProgress(message = paste0('Query for ', ccds, '/', geneS), value = 0.1, {
      
      incProgress(0.2, detail = "(Obtaining continental population)")
        myValue$summary_table <<-
          createMainTable2(geneS, input$depth_of_coverage, per_gene_summary, input$assembly)[,c(1,2,6:10)]

      myValue$assembly = input$assembly
      
      incProgress(0.2, detail = "(Obtaining gene panel tests)")
      myValue$GPT_table <<- createGPT(selectedRow, main_table(), gpt_tP_tG)
      
      incProgress(0.2, detail = "(Creating violin plot)")
      myValue$violon_population <<- createPlot(geneS, main_table()[selectedRow, 2])

      incProgress(0.2, detail = "(Saving details)")
      myValue$gene <<- geneS
      myValue$ccds <<- ccds
      
      setProgress(1)
    })
    showModal(modal_main())

    ##Reset the select_button
    session$sendCustomMessage(type = 'reset_detail_button', message =  "detail_button")

    # reset the "check" on the button in the main table
    session$sendCustomMessage(type = 'resetInputValue', message =  "detail_button")
  })
  
  # create reactive main table
  main_table <- reactive ({
    geneS <- c()
    
    if(length(input$phen) != 0) {
      if (length(input$gene_symbol) != 0) {
        geneS <- c(geneS, input$gene_symbol)
      } else if (length(input$gpt) != 0) {
        geneG <- as.character(unique(gpt_tP_tG[gpt_tP_tG$test_name %in% input$gpt, "gene_symbol"]))
        geneS <- c(geneS, geneG)
      } else if(length(input$phen) != 0) {
        geneP <- unique(gpt_tP_tG$gene_symbol[ 
          gpt_tP_tG$GTR_accession %in% gpt_tP_tG$GTR_accession[ gpt_tP_tG$phenotype_name %in% as.character(input$phen) ]])
        message(length(geneP))
        geneS <- c(geneS, geneP)
        # including HPO
      } #else if(length(input$HPO) != 0) {
       #   geneH <- as.character(unique(ORPHA[ ORPHA$HPO_term_name %in% as.character(input$HPO),"gene_symbol"]))
       #   message(length(geneH))
       #   geneS <- c(geneS, geneH)
       # }
    }
    
    # if statement for HPO
    
    if(length(input$HPO) != 0) {
      if (length(input$gene_symbol) != 0) {
        geneS <- c(geneS, input$gene_symbol)
      } else if(length(input$HPO) != 0) {
           geneH <- as.character(unique(ORPHA[ ORPHA$HPO_term_name %in% as.character(input$HPO),"gene_symbol"]))
           message(length(geneH))
           geneS <- c(geneS, geneH)
      }
    }
    # if (length(input$HPO) != 0) {
    #   geneH <- as.character(unique(ORPHA[ ORPHA$HPO_term_name %in% as.character(input$HPO),"gene_symbol"]))
    #   geneS <- c(geneS, geneH)
    # }
    
    if (length(input$gpt) != 0) {
      if (length(input$phen) != 0 & length(input$HPO) == 0) {
        geneG <- as.character(unique(gpt_tP_tG[gpt_tP_tG$test_name %in% input$gpt, "gene_symbol"]))
        geneS <- c(geneS, geneG)
      }
      if (length(input$phen) == 0 & length(input$HPO) != 0) {
        geneG <- as.character(unique(ORPHA[ORPHA$test_name %in% input$gpt, "gene_symbol"]))
        geneS <- c(geneS, geneG)
      }
    }
    
    if (length(input$gene_symbol) != 0) {
      geneS <- c(geneS, input$gene_symbol)
    }
    
    geneS <- unique(geneS)
    tbl <- createMainTable2(geneS, input$depth_of_coverage, per_gene_summary, input$assembly)
    if(ncol(tbl) > 0) {
      tbl <- tbl[ , c(1:5, 11:13)]
      tbl$Action <- shinyInput(actionButton, nrow(tbl), 'button_', 
                               label = "Detail", 
                               onclick = 'Shiny.onInputChange(\"detail_button\",  this.id)')
    }
    tbl
  })
  
  
  
  observeEvent(input$update, {
    updateNumericInput( session = session, inputId = 'refresh_helper', value = input$refresh_helper + 1 )
  })
  
  # display reactive main table
  output$tableMain <- renderDataTable( {
    #input$update
    message("RH: ", input$refresh_helper)
    t = input$refresh_helper
    if(t > 0 ) {
      isolate({
        if(ncol(main_table()) > 0) {
          formatStyle(datatable(main_table(), escape = FALSE, selection = 'none'), 
                    columns = "Global coverage (mean, %)", target = 'row', backgroundColor = 
                      styleInterval(cuts = c(95, 98), values=c("#FFE4E1", "#FFFFE0", "#F0FFF0")))
        }
      })
    } else {
      data.frame()
    }
  })
  
  # modal window to show tables and plots
  modal_main <- function(failed = FALSE){
    modalDialog(size="l",
      tabsetPanel(type="tabs",
        tabPanel("Per-population summary",
          dataTableOutput('summary_table')),
        tabPanel("Coverage plots",
          fluidRow(
            column(4, align="center", plotOutput("violin_population")),
            column(8, align="center", plotOutput("gnomAD_plot"))
          )
        ),
        tabPanel("Comparison of distribution",
                 fluidPage(
                   sidebarLayout(
                     sidebarPanel(selectInput("KS_1", h3("Select Population 1"),
                                  choices = c("AFR","AMR","EAS","EUR","SAS"),
                                  selected = "AFR"
                                  ),
                                  selectInput("KS_2", h3("Select Population 2"),
                                              choices = c("AFR","AMR","EAS","EUR","SAS"),
                                              selected = "AMR"
                                  )
                                ),
                     mainPanel(
                       plotOutput("KS_plot")
                     )
                   )
                 )
               ),
        tabPanel("Comparison of means",
                 fluidRow(
                   tags$br(),
                   tags$br(),
                 #   column(align = "center",
                    plotOutput("THSD_plot"), height = "100%", width = "100%")
               #   )
        ),
        # including a panel for Tukey results
        
        
        tabPanel("Tests in GTR for the gene",
          dataTableOutput('GPT_table'))
      ),
      easyClose = TRUE
    )
  }

  # table displayed in our modal
  output$summary_table <- renderDataTable({
    myValue$summary_table
  }, escape = FALSE)
  
  output$violin_population <- renderPlot({
    createPlot(myValue$gene, myValue$ccds)
  })
  
  output$GPT_table <- renderDataTable({
    myValue$GPT_table
  }, escape = FALSE)
  
  output$gnomAD_plot <- renderPlot({
    createPlot_gnomAD(myValue$gene)
  })
  
  output$KS_plot <- renderPlot({
    createPlot_KS(myValue$gene, myValue$ccds)
  })
  
  output$THSD_plot <- renderPlot({
    createPlot_THSD(myValue$gene, myValue$ccds)
  })
  
  createPlot <- function(gene, selCCDS) {
    library(dplyr)
    message("[PLOT] Number of CCDS: ", selCCDS, "; Number of genes:", gene, "; Depth:", input$depth_of_coverage )
    message("[boxplot] start -- ", Sys.time())
    selCCDS<-as.character(selCCDS)
    
    idx <- 0 
    dta<-read.fst(sprintf("data/%s_%s.fst", input$assembly, input$depth_of_coverage))
    dta<-melt(filter(dta, ccds_id %in% selCCDS), id.vars = c("ccds_id","gene"))
    dta<-left_join(x=dta, y=kgp.map, by=c("variable" = "ID")) %>% filter(!is.na(pop))
    
    dta$ccds_id <- as.character(dta$ccds_id)
    dta$pop <- as.character(dta$pop)
    dta$variable <- as.character(dta$variable)
    dta$value <- as.numeric(dta$value)

    gAD <- filter(gnomad_exome, ccds_id %in% selCCDS) %>% select(all_of(input$depth_of_coverage)) %>% as.matrix() %>% as.numeric()
    p1 <- ggplot(dta, aes(x=pop, y=value, color=pop)) + 
      theme_bw() + geom_violin()  + 
      facet_wrap(gene~ccds_id) +
      geom_hline(yintercept=gAD, colour="black", show.legend = T) +
      scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
      ylab("Breadth of coverage") + labs(colour = "Population") +
      theme(legend.position="bottom", strip.text.x = element_text(size=12), axis.title=element_text(size=12), axis.title.x = element_blank())
    message("[boxplot]   end -- ", Sys.time())
    p1
  }
  
  # writing new KS plot in modal window
  
  createPlot_KS <- function(gene, selCCDS) {
    message("[PLOT] Number of CCDS: ", selCCDS, "; Number of genes:", gene )
    selCCDS<-as.character(selCCDS)
    
    idx <- 0 
    dta<-read.fst(sprintf("data/%s_%s.fst", input$assembly, input$depth_of_coverage))
    dta<-melt(filter(dta, ccds_id %in% selCCDS), id.vars = c("ccds_id","gene"))
    dta<-left_join(x=dta, y=kgp.map, by=c("variable" = "ID")) %>% filter(!is.na(pop))
    
    dta$ccds_id <- as.character(dta$ccds_id)
    dta$pop <- as.character(dta$pop)
    dta$variable <- as.character(dta$variable)
    dta$value <- as.numeric(dta$value)
    
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    n = 5
    cols = gg_color_hue(n)

    pop_1 <- dta[dta$pop==input$KS_1, "value"]
    colors_1 <-switch(input$KS_1,
                      "AFR" = cols[1],
                      "AMR" = cols[2],
                      "EAS" = cols[3],
                      "EUR" = cols[4],
                      "SAS" = cols[5])
    pop_2 <- dta[dta$pop==input$KS_2, "value"]
    colors_2 <-switch(input$KS_2,
                      "AFR" = cols[1],
                      "AMR" = cols[2],
                      "EAS" = cols[3],
                      "EUR" = cols[4],
                      "SAS" = cols[5])
    group <- c(rep(input$KS_1, length(pop_1)), rep(input$KS_2, length(pop_2)))
    # colors <- c(rep(colors_1, length(pop_1)), rep(colors_2, length(pop_2)))
    
    dat <- data.frame(KSD = c(pop_1,pop_2), group = group)
    
    names(cols) <- levels(dta$group)
    colScale <- scale_colour_manual(name = "group",values = c(colors_1, colors_2))
    
    
    cdf1 <- ecdf(pop_1)
    cdf2 <- ecdf(pop_2)
    
    minMax <- seq(min(pop_1, pop_2), max(pop_1, pop_2), length.out=length(pop_1)) 
    x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )]
    y0 <- cdf1(x0)
    y1 <- cdf2(x0)
    
    test <- ks.test(pop_1, pop_2, alternative="two.sided", exact = NULL)
    plot <- ggplot(dat, aes(x = KSD, group = group, color = group))+
      stat_ecdf(size=1) + colScale +
      theme_bw() +
      theme(legend.position ="top") +
      xlab("Sample") +
      ylab("ECDF") +
      #geom_line(size=1) +
      geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]),
                   linetype = "dashed", color = "black") +
      geom_point(aes(x = x0[1] , y= y0[1]), color="black", size=6) +
      geom_point(aes(x = x0[1] , y= y1[1]), color="black", size=6) +
      theme_update(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0(selCCDS, "\n","(D = ", signif(test$statistic, 3), ", p-value = ",
                     format.pval(test$p.val, digits = 3, eps = 0.001), ")"))
    plot
  }
  
  # writing THSD plot
  
  createPlot_THSD <- function(gene, selCCDS) {
    message("[PLOT] Number of CCDS: ", selCCDS, "; Number of genes:", gene )
    selCCDS<-as.character(selCCDS)
    
    idx <- 0 
    dta<-read.fst(sprintf("data/%s_%s.fst", input$assembly, input$depth_of_coverage))
    dta<-melt(filter(dta, ccds_id %in% selCCDS), id.vars = c("ccds_id","gene"))
    dta<-left_join(x=dta, y=kgp.map, by=c("variable" = "ID")) %>% filter(!is.na(pop))
    
    dta$ccds_id <- as.character(dta$ccds_id)
    dta$pop <- as.character(dta$pop)
    dta$variable <- as.character(dta$variable)
    dta$value <- as.numeric(dta$value)
    
    dta_ANOVA <- aov(dta$value ~ dta$pop)
    dta_THSD <- TukeyHSD(dta_ANOVA, ordered = T)
    
    dta_ANOVA_Population <- dta_THSD$`dta$pop`
    diff <- as.numeric(dta_ANOVA_Population[,1])
    names(diff) <- rownames(dta_ANOVA_Population)
    p_value <- as.numeric(dta_ANOVA_Population[,4])
    names(p_value) <- rownames(dta_ANOVA_Population)

    z <- matrix(
      c(0, diff["AFR-AMR"],diff["AFR-EAS"],diff["AFR-EUR"], diff["AFR-SAS"],
        diff["AMR-AFR"], 0, diff["AMR-EAS"],diff["AMR-EUR"],diff["AMR-SAS"],
        diff["EAS-AFR"], diff["EAS-AMR"], 0, diff["EAS-EUR"], diff["EAS-SAS"],
        diff["EUR-AFR"], diff["EUR-AMR"], diff["EUR-EAS"], 0, diff["EUR-SAS"],
        diff["SAS-AFR"], diff["SAS-AMR"], diff["SAS-EAS"], diff["SAS-EUR"], 0),
      byrow = T,
      ncol = 5,
      nrow = 5,
      dimnames = list(c("AFR", "AMR", "EAS", "EUR", "SAS"),
                      c("AFR", "AMR", "EAS", "EUR", "SAS"))
    )
    
    z2 <- matrix(
      c(1, p_value["AFR-AMR"],p_value["AFR-EAS"],p_value["AFR-EUR"], p_value["AFR-SAS"],
        p_value["AMR-AFR"], 1, p_value["AMR-EAS"],p_value["AMR-EUR"],p_value["AMR-SAS"],
        p_value["EAS-AFR"], p_value["EAS-AMR"], 1, p_value["EAS-EUR"], p_value["EAS-SAS"],
        p_value["EUR-AFR"], p_value["EUR-AMR"], p_value["EUR-EAS"], 1, p_value["EUR-SAS"],
        p_value["SAS-AFR"], p_value["SAS-AMR"], p_value["SAS-EAS"], p_value["SAS-EUR"], 1),
      byrow = T,
      ncol = 5,
      nrow = 5,
      dimnames = list(c("AFR", "AMR", "EAS", "EUR", "SAS"),
                      c("AFR", "AMR", "EAS", "EUR", "SAS"))
    )

    for (i in 1:nrow(z)) {
      for (j in 1:ncol(z)) {
        if (is.na(z[i,j])) {
          z[i,j] <- as.numeric(diff[paste0(rownames(z)[i],"-",colnames(z)[j])])
          if (is.na(z[i,j])) {
            z[i,j] <- as.numeric(diff[paste0(colnames(z)[j],"-",rownames(z)[i])])
          }
          z2[i,j] <- as.numeric(p_value[paste0(rownames(z2)[i],"-",colnames(z2)[j])])
          if (is.na(z2[i,j])) {
            z2[i,j] <- as.numeric(p_value[paste0(colnames(z2)[j],"-",rownames(z2)[i])])
          }
        }
      }
    }
    
    par(mfrow=c(1,2))
    #par(mar=c(5,4,8,2))
    corrplot(z, method="color",is.corr = F, col = brewer.pal(n=5, name="Blues"), type = "lower", title = "Difference of means",mar=c(0,0,1,0))
    corrplot(z2, method="color",is.corr = F, col = brewer.pal(n=5,name="Blues"), type = "lower", title = "p-value", mar = c(0,0,1,0))
    par(mfrow=c(1,1))

  }
  
  createPlot_gnomAD<-function(gene) {
    message("[gnomad] start -- ", Sys.time())
    library(dplyr)
    library(GenomicRanges)
    library(GenomeInfoDb)
    library(GenomeInfoDbData)
    library(RColorBrewer)
    
    message("[gnomad]  load -- ", Sys.time())
    depthL<-c("5x", "10x", "15x", "20x", "25x", "30x", "50x", "100x")
    sample_data = data.frame(
      sample_id = c("005x","010x","015x","020x","025x","030x","050x","100x"),
      bigWig = sprintf("data/%s.bw", depthL),
      scaling_factor = 1,
      stringsAsFactors = F
    )
    track_data = dplyr::mutate(sample_data, track_id = "Coverage", colour_group = sample_id)
    
    selected_tx = dplyr::filter(tx_metadata, gene_name %in% gene, ccds != "") %>% dplyr::select(transcript_id) %>% as.matrix() %>% as.character()

    message("[gnomad]  data -- ", Sys.time())
    pp<-plotCoverage(tx_exons[selected_tx], tx_cdss[selected_tx], tx_metadata, track_data = track_data, rescale_introns = T,
                     fill_palette = brewer.pal(n=8, 'Spectral'), alpha = 0.5, coverage_type = "both", heights = c(0.7, 0.3), return_subplots_list = T)

    message("[gnomad] plot1 -- ", Sys.time())
    plot_cov<-pp$coverage_plot + theme(legend.position = "bottom", legend.title = element_blank()) +
      scale_colour_brewer(palette = "Spectral", labels = c("Over 5x","Over 10x","Over 15x","Over 20x","Over 25x","Over 30x","Over 50x","Over 100x")) +
      scale_fill_brewer(palette = "Spectral", labels = c("Over 5x","Over 10x","Over 15x","Over 20x","Over 25x","Over 30x","Over 50x","Over 100x")) +
      guides(colour=guide_legend(nrow=1))
    ll<-as_ggplot(get_legend(plot_cov))

    message("[gnomad] plot2 -- ", Sys.time())
    if (length(selected_tx) < 5) {
      ((pp$coverage_plot + ylab("Fraction of individuals with coverage \nover X") + theme(axis.text = element_text(size=12))) / pp$tx_structure / ll) +
        plot_layout(heights = c(0.7, 0.25, 0.05)) + plot_annotation(title = sprintf("gnomAD exome coverage for %s", gene))
    } else {
      ((pp$coverage_plot + ylab("Fraction of individuals with coverage \nover X") + theme(axis.text = element_text(size=12))) / pp$tx_structure / ll) +
        plot_layout(heights = c(0.6, 0.35, 0.05)) + plot_annotation(title = sprintf("gnomAD exome coverage for %s", gene))
    }
  }
}

# Run the app ----
shinyApp(ui = ui, server = server)
