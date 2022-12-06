# make interactive plot of contigs in genomes or metagenomes
# last modified 2022-12-06
# v1.2 switch to ggplot and dplyr
# v1.3 add option to highlight contigs with orange diamonds

library(shiny)
library(ggplot2)
library(dplyr)

#coveragefile = "~/project/genomes/ephydatia_muelleri/ASM_HIC_394/Emuelleri_lib001_final_assembly.coverage_gc.tab"
#coveragefile = "~/project/climate_lake_metagenome/climate_lake1_scaffolds.stats.w_genus.tab"
coveragefile = "~/git/genome-reannotations/jbrowse-tracks/sycon/sycon_w_370_600_and_600_700_hisat2.gc_coverage.tab"

coveragedata = read.table(coveragefile, header=TRUE, sep='\t')
#coveragedata = coveragedata[rev(1:nrow(coveragedata)),]

# check if header is correct
if (colnames(coveragedata)[1]!="scaffold"){
  print("ERROR: unrecognized header format in coverage stats file")
  print("headings should be:   scaffold	number	length	coverage	GC	gaps")
}
# contiglengths = coveragedata[["length"]]
# magnituderange = range( log10(contiglengths) )
# pchsize = log10(contiglengths) - magnituderange[1]
# # sizes of three reference points in the legend
# contignames = coveragedata[,1]
# # default is green
# pointcolor = rep("#39bc6744", length(contignames))
# longcontigs = contiglengths > 100000
# massivecontigs = contiglengths > 1000000
# # midsize is blue, must be 100kb
# pointcolor[longcontigs] = "#386edc66"
# # longest contigs are magenta, must be over 1Mb
# pointcolor[massivecontigs] = "#d51ea477"
#point_palette = c("#39bc6744", "#386edc66", "#d51ea477")
point_palette = c(colorRampPalette( c("#81ce9c44", "#81ce9c44", "#386edc77"), alpha=TRUE)(50), 
      colorRampPalette( c("#386edc77", "#860f9b99"), alpha=TRUE)(10), 
      colorRampPalette( c("#860f9b99", "#6e0b62bb"), alpha=TRUE)(40) )

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Contig coverage"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
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
            step = 10
      ),
      sliderInput(inputId = "length",
                  label = "Contig size (log bp)",
                  min = 0,
                  max = 10,
                  value = c(3,8),
                  step = 0.1
      ),
      textInput("userContig", label = "Highlight contigs/scaffolds by name from table, separated by ,", 
                value = ""),
      radioButtons("cov_mode", h4("Coverage axis mode"),
                   choices = c("Linear", "Logarithmic"),
                   selected = "Linear"
                   ),
      h4("Render current graph as PDF"),
      downloadButton("printpdf", label = "Print")
    ),
    # Main panel for displaying outputs ----
    mainPanel( 
      h3("Each point is a contig. Click-and-drag to display stats"),
      #strong( paste("Using", basename(coveragefile) ) ),
      plotOutput(outputId = "distPlot",
                 height="600px",
                 click = "plot_click",
                 brush = brushOpts(id = "plot_brush")
                 ),
      tableOutput("info")
      )
   )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  make_covplot = reactive({
    # axis_mode = ""
    # coverage_range = input$cov
    # if (input$cov_mode == 2){
    #   axis_mode = "x"
    #   if ( input$cov[1]==0 ){
    #     coverage_range = c(1, input$cov[2] )
    #   }
    # }

    # ggplot tries to resize and recolor based on the filtered data
    # setting the variables here makes the values objective between datasets
    # allowing for comparison
    fd = mutate(coveragedata, #[rev(1:nrow(coveragedata)),],
                lencolor = point_palette[round(log10(coveragedata$length)*10)],
                lensize = log(coveragedata$length, base=8)
                ) %>%
      filter( log10(length) >= input$length[1] &
              log10(length) <= input$length[2]
                )  %>%
      arrange(length) # put biggest points last, so they appear on top
    
    uf = filter(coveragedata, scaffold %in% as.list(scan(text=trimws( input$userContig ), what='', sep=',') ) )
    
    all_contig_total = sum(coveragedata[["length"]])
    sub_table = brushedPoints(coveragedata, input$plot_brush, xvar = "coverage", yvar = "GC")
    sub_size_total = sum(sub_table[["length"]])

    gg = ggplot(fd, aes( x=coverage, y=GC , color=lencolor ) ) +
      theme_minimal() +
      theme(legend.position="none",
            axis.text=element_text(size=16),
            axis.title=element_text(size=18),
            plot.title = element_text(size=20) ) +
      scale_y_continuous(limits=input$gc) +
      scale_x_continuous(expand = expansion(mult = c(0.005, 0.03) ),
                         limits=c( ifelse(input$cov_mode=="Linear",input$cov[1],1), input$cov[2]),
                         trans = ifelse(input$cov_mode=="Linear","identity","log10") ) +
      scale_colour_identity() +
      #scale_colour_stepsn(colours = point_palette, breaks = c(0, 100000, 1000000) ) +
      labs(x="Coverage", y="GC%", title=paste("Using", basename(coveragefile) ) ) +
      geom_point( size=fd$lensize ) +
      annotate(geom="text", x=min(input$cov), y=max(input$gc), size=5, hjust=0,
               label=paste(dim(coveragedata)[1], "contigs", round(all_contig_total/1000000,digits=1),"Mb total") ) + 
      annotate(geom="text", x=max(input$cov), y=min(input$gc), size=5, hjust=1,
               label=paste(dim(sub_table)[1], "contigs", round(sub_size_total/1000000,digits=1),"Mb selected") ) +
      annotate(geom="text", x=max(input$cov), y=max(input$gc), size=5, hjust=1,
               label=paste("Showing", nrow(fd), "contigs from", round(10^input$length[1]), "bp to", round(10^input$length[2]), "bp") )
    # add optional highlighting points
    if (TRUE){
      # draw black points for all bounded box contigs
      gg = gg + geom_point( data=sub_table, col="#000000" ) 
      # show user picked scaffold, if any
      gg = gg + geom_point( data = uf, aes(x=coverage, y=GC), color ="#f8520a", size=7, stroke=3, shape=5, alpha=0.8)
    }
    
    gg
    
  # leftover code from original script
    # datarange = input$length[1] <= log10(contiglengths) & log10(contiglengths) <= input$length[2]
    # par(mar=c(4.5,4.5,3,1))
    # # reverse row order, since it usually expects biggest contigs first
    # # this ends up with largest contigs plotted last, meaning top layer
    # plot( rev(coveragedata[["coverage"]][datarange]), rev(coveragedata[["GC"]][datarange]), type='p',
    #      xlim=coverage_range, ylim=input$gc, log=axis_mode,
    #      xlab="Mean coverage of mapped reads", ylab="GC%",
    #      pch=16, frame.plot=FALSE, col=rev(pointcolor[datarange]), cex.axis=1.5, cex=rev(pchsize[datarange]), main="", cex.lab=1.4)
    # # display overview stats on top left of graph
    # text(min(input$cov), max(input$gc), paste(dim(coveragedata)[1], "contigs", round(all_contig_total/1000000,digits=1),"Mb total"), pos=4, cex=1.2)
    # # display selected contigs on bottom right of graph
    # text(max(input$cov), min(input$gc), paste(dim(sub_table)[1], "contigs", round(sub_size_total/1000000,digits=1),"Mb selected"), pos=2, cex=1.2)
    # legend(covmax,gcmax, legend=legendlabels, pch=16, col=c("#39bc6799","#386edc99","#d51ea499"), pt.cex=legendpch, cex=1.1, title="Contig size (bp)", xjust=1)
    # # text(covmax,23,paste(round(totalsize/1000000,digits=1),"Mb"), cex=1.2, pos=2)
    # # text(covmax,20,paste(length(contignames),"total contigs"), cex=1.2, pos=2)
  })
  
  output$distPlot <- renderPlot({
    make_covplot()
  })
  
  output$info <- renderTable({
    brushedPoints(coveragedata, input$plot_brush, xvar = "coverage", yvar = "GC")
    # previous version, where points were clicked
    #nearPoints(coveragedata, input$plot_click, xvar = "coverage", yvar = "GC")
  })
  
  output$printpdf <- downloadHandler(
    filename = function() {"plot.pdf"},
    content = function(filename){
      gg_covplot = make_covplot()
      ggsave(filename, gg_covplot, device="pdf", width=8, height=6)
    }
  )
  # output$printpdf <- downloadHandler(
  #   filename = function() {"plot.pdf"},
  #   content = function(file) {
  #     axis_mode = ""
  #     coverage_range = input$cov
  #     if (input$cov_mode == 2){
  #       axis_mode = "x"
  #       if ( input$cov[1]==0 ){
  #         coverage_range = c(1, input$cov[2] )
  #       }
  #     }
  #     datarange = input$length[1] <= log10(contiglengths) & log10(contiglengths) <= input$length[2]
  #     pdf_title = paste0(input$cov[1],"-",input$cov[2],"cov_", input$gc[1], "-",input$gc[2], "gc", ".pdf")
  #     pdf(file, width=7, height=6, title=pdf_title)
  #     par(mar=c(4.5,4.5,3,1))
  #     plot( rev(coveragedata[["coverage"]][datarange]), rev(coveragedata[["GC"]][datarange]), type='p', 
  #          xlim=coverage_range, ylim=input$gc, log=axis_mode,
  #          xlab="Mean coverage of mapped reads", ylab="GC%", 
  #          pch=16, frame.plot=FALSE, col=rev(pointcolor[datarange]), cex.axis=1.5, cex=rev(pchsize[datarange]), main="", cex.lab=1.4)
  #     dev.off()
  #     }
  #   )
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)




