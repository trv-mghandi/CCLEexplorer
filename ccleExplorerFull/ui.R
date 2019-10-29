## codepSearch
#if(TRUE){

topdepfluidRow2 <- fluidRow(
  tags$script("
    Shiny.addCustomMessageHandler('doubleClick', function(value) {
    Shiny.setInputValue('doubleClick', value);
    });
  ")
    ,
    box(
    width = 6, status = "info", solidHeader = TRUE,
    title = "correlated features network"
    #,
    #bubblesOutput("packagePlot", width = "100%", height = 600)
    #plotlyOutput("plot", height = 900)#, width = "100%")#, height = 800)
    ,fluidRow(
      column(3,
             selectInput(inputId = "query.feature.type",#inputId = "query.gene",
                         label = "dataset",
                         choices = f.types.nam, 
                         selected = f.types.nam[grep('AVANA', toupper(f.types.nam))[1]])
          #              actionButton("refresh.plot.network", label = "Refresh plots") 
      ),
      
      column(3,
             selectInput(inputId = "query.feature.name",
                         label = "feature",
                         choices = colnames(F)[F.type=="Avana"],#f.fnams.in.type, 
                         selected = 'Avana:ESR1'
                         #              actionButton("refresh.plot.network", label = "Refresh plots") 
             )),

      column(3,
               selectInput(inputId = "query.nmax",
                           label = "# top features",
                           choices = nmax.options, 
                           selected = 20) ) 
      # ,column(3,
      #           # Input: Numeric entry for number of obs to view ----
      #           selectInput(inputId = "query.dmax",
      #                       label = "max dist",
      #                       choices = dmax.options, 
      #                       selected = 1)
      # ) 
      ,      column(3,
                selectInput(inputId = "query.emax",
                            label = "max edges shown",
                            choices = emax.options, 
                            selected = 20)
      ))
    # , 
    # fluidRow(
    #   column(3,
    #          selectInput(inputId = "query.feature.name",
    #                      label = "feature",
    #                      choices = colnames(F),#f.fnams.in.type, 
    #                      selected = 'Avana:ESR1'
    #                      #              actionButton("refresh.plot.network", label = "Refresh plots") 
    #          )))
    # 
    ,visNetworkOutput("networkTopCodep2",height = 600)
    ,checkboxGroupInput("Filt.types", "Feature types to show:",
                        g.cur.f.types.nam, inline = TRUE, selected = g.cur.f.types.nam.selected), 
    
    fluidRow(
      column(2,actionButton("check.all", label = "check all")),
      column(2,actionButton("uncheck.all", label = "uncheck all"))
    ) 
    
    ,checkboxGroupInput("add.filts", "Additional filters:", c('related.genes', 'remove.redundancies'), inline = TRUE
    )
    
    ,verbatimTextOutput("shiny_return")
    #,tableOutput("data")
    
  )
  , box(
    width = 6, status = "info", solidHeader = TRUE,
    title = "scatter plot",
    
    fluidRow(
      column(3,
             selectInput(inputId = "x.dataset",
                         label = "x axis dataset",
                         choices = f.types.id)#g$dataSets)
             #              actionButton("refresh.plot.network", label = "Refresh plots") 
      ),
      column(3,
             selectInput(inputId = "x.gene",
                         label = "x axis feature",
                         choices = f.types.id)#g$dataSets)
      ),
      column(3,
             selectInput(inputId = "y.dataset",
                         label = "y axis dataset",
                         choices = f.types.id)#g$dataSets)
      ),
      column(3,
             selectInput(inputId = "y.gene",
                         label = "y axis feature",
                         choices = f.types.id)#g$dataSets)
      )
      
    ),
    
    plotlyOutput("plot", height = 600)#, width = "100%")#, height = 800)
    ,
    fluidRow(
      column(2,
             selectInput(inputId = "z.dataset",
                         label = "z axis dataset",
                         choices = f.types.id)#g$dataSets)
      ),
      column(2,
             selectInput(inputId = "z.gene",
                         label = "z axis feature",
                         choices = f.types.id)#g$dataSets)
      ),
      column(3,selectInput(inputId = "color.by",
                           label = "color by",
                           choices = c('',g$color_type_family, 'Custom'), selected = 'type_refined')), 
      
      column(2,
             actionButton("refresh.plot", label = "Refresh plot")),
      
      column(1,actionButton("toggle.2D3D", label = "2D/3D")),
      
      column(2,actionButton("toggle.all", label = "toggle.all"))
      
    )  
  ) 
  
  
) 


dashboardPage(
  dashboardHeader(title = paste0("CCLE Explorer",ifelse(DEMO_mode,' DEMO',''))), 
  dashboardSidebar(collapsed = TRUE,
                   #  sliderInput("rateThreshold", "Warn when rate exceeds",
                   #              min = 0, max = 50, value = 3, step = 0.1
                   #  ),
                   sidebarMenu(
                     #menuItem("Top co-dependencies", tabName = "codep_top"),
                     menuItem("Search co-dependencies", tabName = "codep_one")#,
                   )
  ),
  dashboardBody(
    tabItems(
      
      tabItem("codep_one",
              topdepfluidRow2
      )#,
      
      #tabItem("codep_one",
      #       box(
      #        width = 10, status = "info", solidHeader = TRUE,
      #       title = "cdep_one"#,
      #bubblesOutput("packagePlot", width = "100%", height = 600)
      #plotlyOutput("plot", height = 900)#, width = "100%")#, height = 800)
      #    )
      #)
    )
  )
)

#}
