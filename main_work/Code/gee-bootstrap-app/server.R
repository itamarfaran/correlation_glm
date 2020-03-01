library(data.table)
library(DT)
library(plotly)

data <- fread('gee_data.csv')
data[,`:=`(
  p_n_ratio = p/n,
  nh = sick_obs_percentage*n,
  nd = (1 - sick_obs_percentage)*n
)]
data_plt <- copy(data)

data[,`:=`(
  actual_sd = round(actual_sd, 3),
  gee_sd = round(gee_sd, 3),
  gee_old_sd = round(gee_old_sd, 3),
  mle_sd = round(mle_sd, 3)
)]
data_plt[,`:=`(
  none_ = factor(1),
  p_ = factor(p),
  n_ = factor(n),
  p_n_ratio_ = p_n_ratio,
  nh_ = factor(ceiling(nh)),
  nd_ = factor(ceiling(nd)),
  sick_obs_percentage_ = factor(sick_obs_percentage),
  Tlength_ = factor(Tlength)
  )]

function(input, output){
  output$points <- renderPlot({
    data_plt[,ratio := get(input$y_numerator)]
    if(input$y_denumerator != '1') data_plt[,ratio := ratio/get(input$y_denumerator)]
    out <- ggplot(data_plt, aes_string(
      x = input$x,
      y = 'ratio',
      col = paste0(input$color, '_')
    )) +
      geom_point(position = (if(input$jitter) 'jitter' else 'identity')) +
      geom_hline(yintercept = 0, size = 1) +
      labs(
        x = sidebar_options_reverse[input$x],
        y = paste0(sidebar_options_reverse[input$y_numerator], ' / ' , sidebar_options_reverse[input$y_denumerator]),
        col = sidebar_options_reverse[input$color]
      )
    if(!(input$aggFun == 'none')){
      out <- out + switch(
        input$aggFun,
        'mean' = stat_summary(fun.y = mean, geom = "line", size = 1.2, linetype = 1),
        'median' = stat_summary(fun.y = median, geom = "line", size = 1.2, linetype = 1),
        'smooth' = geom_smooth(se = FALSE),
        'lm' = geom_smooth(se = FALSE, method = 'lm'))
    }
    if(input$vline > 0) out <- out + geom_hline(yintercept = as.numeric(input$vline), size = 1, col = 'darkgrey', linetype = 2)
    if(input$abline) out <- out + geom_abline(slope = 1, intercept = 0, size = 1, col = 'darkgrey', linetype = 2)
    out
    })
  output$boxplot <- renderPlot({
    data_plt[,ratio := get(input$y_numerator)]
    if(input$y_denumerator != '1') data_plt[,ratio := ratio/get(input$y_denumerator)]
    out <- ggplot(data_plt, aes_string(
      x = input$x,
      group = input$x,
      y = 'ratio'
    )) +
      geom_boxplot(fill = 'lightblue') +
      geom_hline(yintercept = 0, size = 1) +
      labs(
        x = sidebar_options_reverse[input$x],
        y = sidebar_options_reverse[input$y_numerator]
        )
    if(input$y_denumerator != '1') out <- out + geom_hline(yintercept = 1, size = 1, col = 'darkgrey', linetype = 2)
    out
  })
  output$data <- DT::renderDataTable(data)
}