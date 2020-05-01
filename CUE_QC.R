#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
#library(datasets)

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Quality Control for CUE"),
  tags$h5("This app allows one to select his/her own QC criteira for whole blood or placenta samples. The default QC criteira for whole blood samples are MSE<0.0025 and Accruracy > 95%. The defalut QC criteria for placenta samples are MSE<0.01 and Accruracy > 90%."),
  
  fluidRow(column(3, 
                  radioButtons("radio", 
                               h3("Tissues"), 
                               choices = list("1. Whole Blood (PTSD)" = "Whole Blood (PTSD)", 
                                              "2. Placenta (ELGAN)" = "Placenta (ELGAN)"),
                               selected = "Whole Blood (PTSD)")),
  ),
  
  textOutput("tissueText"),
  
  htmlOutput("QCText"),
  tags$br(),
  
  htmlOutput("SummaryText"),
  
  tags$h3("Click the button below to download the QC+ CpG list!"),
  downloadButton("downloadData", "Download"),
  
  tags$br(),
  tags$br(),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("MSE_Bar",
                  "MSE_Sep",
                  min = 0,
                  max = 0.01,
                  value = 0.0025)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("MSEPlot")
    )
  ),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("Accuracy_Bar",
                  "Accuracy_Sep",
                  min = 0.9,
                  max = 1,
                  value = 0.95)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("accPlot")
    )
  ),
  
  
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$tissueText <- renderText({
    paste("You have selected dataset ", input$radio,".",sep="")
  })
  
  
  output$MSEPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    if (input$radio=="Whole Blood (PTSD)"){
      load("Data/PTSD_CpG_Best_method_list.RData")
      load("Data/PTSD_CUE_Evaluations.RData")
    } else{
      load("Data/ELGAN_CpG_Best_method_list.RData")
      load("Data/ELGAN_CUE_Evaluations.RData")
    }
    
    x    <- mse.min
    MSE_Bar <- seq(0, 0.01, length.out = 20000 + 1)
    #bins <- seq(min(x), max(x), length.out = input$MSE_Bar + 1)
    
    # draw the histogram with the specified number of bins
    #hist(x, breaks = bins, col = 'darkgray', border = 'white')
    plot(sort(x),type="l",#col='cyan',
         xlab = "Number of Sites",
         ylab = "MSE",
         main = "MSE Threshold Selection")
    abline(h=input$MSE_Bar,col="red")
    
  })
  
  output$accPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    if (input$radio=="Whole Blood (PTSD)"){
      load("Data/PTSD_CpG_Best_method_list.RData")
      load("Data/PTSD_CUE_Evaluations.RData")
    } else{
      load("Data/ELGAN_CpG_Best_method_list.RData")
      load("Data/ELGAN_CUE_Evaluations.RData")
    }
    x    <- acc.max
    MSE_Bar <- seq(0.9, 1, length.out = 100 + 1)
    #bins <- seq(min(x), max(x), length.out = input$MSE_Bar + 1)
    
    # draw the histogram with the specified number of bins
    #hist(x, breaks = bins, col = 'darkgray', border = 'white')
    plot(sort(x, decreasing = TRUE),type="l",#col='cyan',
         xlab = "Number of Sites",
         ylab = "Accuracy",
         main = "Accuracy Threshold Selection")
    abline(h=input$Accuracy_Bar,col="red")
    
  })
  
  output$QCText <- renderUI({
    str1 <- paste("You have selected the MSE threshold: ", input$MSE_Bar)
    str2 <- paste("You have selected the accuracy threshold: ", input$Accuracy_Bar)
    HTML(paste(str1, str2, sep = '<br/>')) 
  })
  
  output$SummaryText <- renderUI({
    if (input$radio=="Whole Blood (PTSD)"){
      load("Data/PTSD_CpG_Best_method_list.RData")
      load("Data/PTSD_CUE_Evaluations.RData")
    } else{
      load("Data/ELGAN_CpG_Best_method_list.RData")
      load("Data/ELGAN_CUE_Evaluations.RData")
    }
    x <- acc.max
    y <- mse.min
    str1 <- paste("Number of sites below MSE threshold: ",  sum(y < input$MSE_Bar), " (",round(100*sum(y < input$MSE_Bar)/length(y),2),"%)",sep="")
    str2 <- paste("Number of sites above the accuracy threshold: ",  sum(x > input$Accuracy_Bar), " (",round(100*sum(x > input$Accuracy_Bar)/length(y),2),"%)",sep="")
    str3 <- paste("Number of sites passed QC: ",  sum((x > input$Accuracy_Bar)+( y < input$MSE_Bar)==2 ), " (",round(100*sum((x > input$Accuracy_Bar)+( y < input$MSE_Bar)==2 )/length(y),2),"%)",sep="")
    HTML(paste(str1, str2,str3, sep = '<br/>'))
    #HTML(paste("MSE threshold: ", input$MSE_Bar, "; Accuracy threshold: ", input$Accuracy_Bar,".",sep=""))
    
  })
  
  # Reactive value for selected dataset ----
  # datasetInput <- reactive({
  #     switch(input$dataset,
  #            "rock" = rock,
  #            "pressure" = pressure,
  #            "cars" = cars)
  # })
  
  # Table of selected dataset ----
  output$table <- renderTable(which((x > input$Accuracy_Bar)+( y < input$MSE_Bar)==2 ))
  #     renderTable({
  #     datasetInput()
  #     
  # })
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    
    filename = function() {
      paste(input$radio, "QC_probe_list.csv", sep = "")
    },
    content = function(file) {
      if (input$radio=="Whole Blood (PTSD)"){
        load("Data/PTSD_CpG_Best_method_list.RData")
        load("Data/PTSD_CUE_Evaluations.RData")
      } else{
        load("Data/ELGAN_CpG_Best_method_list.RData")
        load("Data/ELGAN_CUE_Evaluations.RData")
      }
      x <- acc.max
      y <- mse.min
      z <- names(mse.min)
      #write.csv(which((x > input$Accuracy_Bar)+( y < input$MSE_Bar)==2 ), file, row.names = FALSE)
      write.csv(z[which((x > input$Accuracy_Bar)+( y < input$MSE_Bar)==2) ], file, row.names = FALSE)
    }
  )
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
