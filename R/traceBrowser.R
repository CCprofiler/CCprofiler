#' traceBrowser for SECprofiler
#' @description basic data browser allowing to display select subsets of chromatographic elution
#' profiles, e.g. SECtraces
#' @param Traces traces object, peptide or protein level
#' @param calibration_functions calibration object created by function calibration, providing functions to map
#' protein_mw to fraction and back, allowing to indicate expected monomeric elution fraction of selected proteins
#' @import data.table
#' @import plotly
#' @import shiny
#' @import reshape2
#' @import ggplot2
#' @export

traceBrowser <- function(Traces, calibration_functions = NULL)
  {
  toLongFormat <- function(traces.dt) {
  traces.dt.long <-
    melt(traces.dt, id.var='id', variable.name='fraction',
         value.name='intensity', variable.factor=FALSE)
  traces.dt.long[, fraction := as.numeric(fraction)]
  setkey(traces.dt.long,id)
  data.table(traces.dt.long)
  traces.dt.long
}
  # Prepare data for input fields
  ################################
  # get annotations and assign vectors as fvalue input
  annotations <- names(Traces$trace_annotation)
  for (i in seq_along(annotations)){
    assign(annotations[i], Traces$trace_annotation[[i]])
  }
  
  # build and run Shiny App to visualize subset traces
  shinyApp(
    ui = pageWithSidebar(
      
      # Application title
      headerPanel("SECprofiler Trace Browser"),
      
      # Sidebar with a slider input for number of observations
      sidebarPanel(
        selectInput("fcolumn", "Choose Annotation Column for filtering", annotations, selected = "id"),
        
        uiOutput("fcolumnvalues"), #The possible choices of this field are calculated on the server side and passed over by uiOutput
        
        textInput("legend.position", label = "Plot legend position (or \"none\")", value = "right"),
        
        p("Select annotation column and value to display matching traces (i.e. elution profiles).
          dotted vertical lines give theoretical monomeric elution fractions of individual proteins")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotlyOutput("plot"),
        tableOutput("table")
      )
    ), 
    server = function(input, output) {
      
      ## Generate Reactive Filter Value Field for UI, depending on filter column chosen
      output$fcolumnvalues <- renderUI({
        values <- sort(unique(get(input$fcolumn)))
        # values <- values[nchar(values)>0]
        selectizeInput("fvalue", "Choose Filter Value(s)", values, multiple = TRUE, options = list(maxOptions = 60000), selected = values[1])
      })
      
      ## generate selected protein SEC traces plot
      output$plot <- renderPlotly({
        
        # Collect indices of the proteins matching the selection 
        indices <- numeric()
        for (i in 1:length(input$fvalue)){
          column_index <- which(annotations == input$fcolumn)
          indices = c(indices, grep(input$fvalue[i], get(input$fcolumn)))
        }
        
        # subset the data for these, including Mass column
        Traces.s <- subset(Traces, trace_ids = Traces$traces$id[indices])
        
        # Plot
        traces.long <- toLongFormat(Traces.s$traces)
        traces.long <- merge(traces.long, Traces.s$trace_annotation, by = "id", all.y = FALSE, all.x = TRUE)
        p <- ggplot(traces.long) +
          ggtitle(paste(Traces.s$trace_type, 'traces selected')) +
          xlab('fraction') + ylab('intensity') + 
          theme_bw() +
          theme(legend.position=input$legend.position) 
        if (is.null(calibration_functions)){
          q <- p + geom_line(aes(x=fraction, y=intensity, color=id))
          ggplotly(q)
        } else {
          lx.frc <- seq(10,(ncol(Traces.s$traces)-1),10)
          lx <- paste( lx.frc , round(calibration_functions$SECfractionToMW(lx.frc), 1) , sep = '\n' )
          
          q <- p + geom_line(aes(x=fraction, y=intensity, color=id)) +
            geom_vline(aes(xintercept = calibration_functions$MWtoSECfraction(protein_mw), color = id , show.legend=FALSE)) +
            scale_x_continuous(breaks = lx.frc , labels = lx) + xlab("fraction\n(apparent MW [kDa])")
          ggplotly(q)
        }
     })
       
      output$table <- renderTable({
       # Collect indices of the proteins matching the selection 
        indices <- numeric()
        for (i in 1:length(input$fvalue)){
          column_index <- which(annotations == input$fcolumn)
          indices = c(indices, grep(input$fvalue[i], get(input$fcolumn)))
        }
        
        Traces.s <- subset(Traces, trace_ids = Traces$traces$id[indices])
        # and subset the data for these, including Mass column
        Traces.s$trace_annotation
      })
    }
  )
}
