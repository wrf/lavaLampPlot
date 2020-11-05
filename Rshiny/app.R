#

library(shiny)

#coveragefile = "~/genomes/ephydatia_muelleri/ASM_HIC_394/Emuelleri_lib002_final_assembly.coverage_gc.tab"
coveragefile = "~/genomes/ephydatia_muelleri/ASM_HIC_394/Emuelleri_lib002_final_assembly.coverage_gc.tab"

coveragedata = read.table(coveragefile, header=TRUE, sep='\t')
contiglengths = coveragedata[["length"]]
magnituderange = range( log10(contiglengths) )
pchsize = log10(contiglengths) - magnituderange[1]
# sizes of three reference points in the legend
contignames = coveragedata[,1]

# default is green
pointcolor = rep("#39bc6744", length(contignames))
longcontigs = contiglengths > 100000
massivecontigs = contiglengths > 1000000
# midsize is blue
pointcolor[longcontigs] = "#386edc66"
# longest contigs are magenta
pointcolor[massivecontigs] = "#d51ea477"


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Contig coverage"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
    #  helpText("Each point is a contig"),
    
      sliderInput(inputId = "gc",
            label = "GC %",
            min = 0,
            max = 100,
            value = c(20,80)
      ),
      sliderInput(inputId = "cov",
            label = "Coverage",
            min = 0,
            max = 2500,
            value = c(0,1000),
            step = 100
      ),
      sliderInput(inputId = "length",
                  label = "Contig size (log bp)",
                  min = 0,
                  max = 9,
                  value = c(1,8),
                  step = 1
      )
    ),
    # Main panel for displaying outputs ----
    mainPanel( 
      h3("Each point is a contig. Click on a point to display stats"),
      strong( paste("Using", basename(coveragefile) ) ),
      plotOutput(outputId = "distPlot",
                 height="600px",
                 click = "plot_click"),
      tableOutput("info")
      )
   )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  output$distPlot <- renderPlot({
    datarange = input$length[1] <= log10(contiglengths) & log10(contiglengths) <= input$length[2]
    par(mar=c(4.5,4.5,3,1))
    plot(coveragedata[["coverage"]][datarange], coveragedata[["GC"]][datarange], type='p', 
         xlim=input$cov, ylim=input$gc, 
         xlab="Mean coverage of mapped reads", ylab="GC%", 
         pch=16, frame.plot=FALSE, col=pointcolor[datarange], cex.axis=1.5, cex=pchsize[datarange], main="", cex.lab=1.4)
   # legend(covmax,gcmax, legend=legendlabels, pch=16, col=c("#39bc6799","#386edc99","#d51ea499"), pt.cex=legendpch, cex=1.1, title="Contig size (bp)", xjust=1)
  #  text(covmax,23,paste(round(totalsize/1000000,digits=1),"Mb"), cex=1.2, pos=2)
  #  text(covmax,20,paste(length(contignames),"total contigs"), cex=1.2, pos=2)
  })
  output$info <- renderTable({
    # With base graphics, need to tell it what the x and y variables are.
    nearPoints(coveragedata, input$plot_click, xvar = "coverage", yvar = "GC")
    # nearPoints() also works with hover and dblclick events
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)













