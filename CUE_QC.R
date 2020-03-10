#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Quality Control for CUE"),

    fluidRow(column(3, 
                    radioButtons("radio", 
                               h3("Tissues"), 
                               choices = list("Whole Blood (PTSD)" = 1, 
                                              "Placenta (ELGAN)" = 2),
                               selected = 1)),
        
        
        
    ),
    textOutput("tissueText"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("RMSE_Bar",
                        "RMSE_Sep",
                        min = 0,
                        max = 0.1,
                        value = 0.0025)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("rmsePlot")
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
    
    textOutput("QCText")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$tissueText <- renderText({
        paste("You have selected", input$radio)
    })
    
    output$rmsePlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        if (input$radio==1){
            load("/Users/GangLi/Dropbox/DNA Methylation Imputation/Code/CUE/PTSD_CpG_Best_method_list.RData")
            load("/Users/GangLi/Dropbox/DNA Methylation Imputation/Code/CUE/PTSD_CUE_Evaluations.RData")
        } else{
            load("/Users/GangLi/Dropbox/DNA Methylation Imputation/Code/CUE/ELGAN_CpG_Best_method_list.RData")
            load("/Users/GangLi/Dropbox/DNA Methylation Imputation/Code/CUE/ELGAN_CUE_Evaluations.RData")
        }
        
        x    <- mse.min
        RMSE_Bar <- seq(0, 0.1, length.out = 200 + 1)
        #bins <- seq(min(x), max(x), length.out = input$RMSE_Bar + 1)

        # draw the histogram with the specified number of bins
        #hist(x, breaks = bins, col = 'darkgray', border = 'white')
        plot(sort(x),type="l",#col='cyan',
             xlab = "Number of Sites",
             ylab = "RMSE",
             main = "RMSE Threshold Selection")
        abline(h=input$RMSE_Bar,col="red")
        
    })
    
    output$accPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        if (input$radio==1){
            load("/Users/GangLi/Dropbox/DNA Methylation Imputation/Code/CUE/PTSD_CpG_Best_method_list.RData")
            load("/Users/GangLi/Dropbox/DNA Methylation Imputation/Code/CUE/PTSD_CUE_Evaluations.RData")
        } else{
            load("/Users/GangLi/Dropbox/DNA Methylation Imputation/Code/CUE/ELGAN_CpG_Best_method_list.RData")
            load("/Users/GangLi/Dropbox/DNA Methylation Imputation/Code/CUE/ELGAN_CUE_Evaluations.RData")
        }
        x    <- acc.max
        RMSE_Bar <- seq(0.9, 1, length.out = 100 + 1)
        #bins <- seq(min(x), max(x), length.out = input$RMSE_Bar + 1)
        
        # draw the histogram with the specified number of bins
        #hist(x, breaks = bins, col = 'darkgray', border = 'white')
        plot(sort(x, decreasing = TRUE),type="l",#col='cyan',
             xlab = "Number of Sites",
             ylab = "Accuracy",
             main = "Accuracy Threshold Selection")
        abline(h=input$Accuracy_Bar,col="red")
        
    })
    
    output$QCText <- renderText({
        paste("RMSE threshold:", input$RMSE_Bar, "; \n Accuracy threshold:", input$Accuracy_Bar,".")
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
