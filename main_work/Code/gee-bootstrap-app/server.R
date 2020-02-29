library(data.table)
library(DT)
library(plotly)
data <- fread('gee_data.csv')
data[,none := factor(1)]

function(input, output){
  output$plt <- renderPlot({
    out <- ggplot(data, aes_string(
      x = input$x,
      y = input$y,
      col = factor(input$color),
      shape = factor(input$shape)
    )) +
      geom_point() +
      geom_hline(yintercept = 0, size = 1) +
      stat_summary(fun.y = mean, geom = "line", size = 1.2, linetype = 1) +
      labs(
        x = sidebar_options_reverse[input$x],
        y = sidebar_options_reverse[input$y],
        col = sidebar_options_reverse[input$color],
        shape = sidebar_options_reverse[input$color]
      )
    if(input$y == 'actual_est_ratio') out <- out + geom_hline(yintercept = 1, size = 1, col = 'darkgrey', linetype = 2)
    out
    })
  output$data <- DT::renderDataTable(data)
}