library(data.table)
library(DT)
library(plotly)
library(shiny)

data <- fread('gee_data.csv')
data[,`:=`(
  none = factor(1),
  m =  0.5*p*(p - 1),
  p_n_ratio = round(p/n, 1),
  m_n_ratio = round(0.5*p*(p - 1)/n, 1),
  nh = 5*round(sick_obs_percentage*n / 5),
  nd = 5*round((1 - sick_obs_percentage)*n / 5)
)]
data_plt <- copy(data)
cols <- colnames(data_plt)
data_plt[,(paste0(cols, '_')) := lapply(.SD, factor), .SDcols = cols]

data[,`:=`(
  actual_sd = round(actual_sd, 3),
  gee_sd = round(gee_sd, 3),
  gee_old_sd = round(gee_old_sd, 3),
  mle_sd = round(mle_sd, 3)
)]


function(input, output){
  output$graph <- renderPlot({
    data_plt[,ratio := get(input$y_numerator)]
    if(input$y_denumerator != '1') data_plt[,ratio := ratio/get(input$y_denumerator)]
    if(input$y_numerator == '1') data_plt[,ratio := 1/ratio]
    if('boxplot' %in% input$plot_checkboxs){
      out <- ggplot(data_plt, aes_string(
        x = input$x,
        group = input$x,
        y = 'ratio'
      )) +
        geom_boxplot(fill = 'lightblue') +
        geom_hline(yintercept = 0, size = 1) +
        labs(
          x = sidebar_options_reverse[input$x],
          y = paste0(sidebar_options_reverse[input$y_numerator], ' / ' , sidebar_options_reverse[input$y_denumerator])
        )
      if(input$y_denumerator != '1') out <- out + geom_hline(yintercept = 1, size = 1, col = 'darkgrey', linetype = 2)
    } else {
      out <- ggplot(data_plt, aes_string(
        x = input$x,
        y = 'ratio',
        col = paste0(input$color, ifelse('color_as_factor' %in% input$plot_checkboxs, '_', ''))
      )) +
        geom_point(position = (if('jitter' %in% input$plot_checkboxs) 'jitter' else 'identity')) +
        geom_hline(yintercept = 0, size = 1) +
        # scales_x[[input$x_scale]] + scales_y[[input$y_scale]] + 
        labs(
          x = sidebar_options_reverse[input$x],
          y = paste0(sidebar_options_reverse[input$y_numerator], ' / ' , sidebar_options_reverse[input$y_denumerator]),
          col = sidebar_options_reverse[input$color]
        )
      if(!(input$aggFun == 'none')){
        x <- 'x'
        # formula <- paste0('y ~ ', ifelse(input$intercept, '', '0 + '), 'poly(x, ', input$polynom, ')')
        formula <- paste0('y ~ poly(x, ', input$polynom, ')')
        formula <- as.formula(formula)
        
        out <- out + switch(
          input$aggFun,
          'mean' = stat_summary(fun.y = mean, geom = "line", size = 1.2, linetype = 1),
          'median' = stat_summary(fun.y = median, geom = "line", size = 1.2, linetype = 1),
          'smooth' = geom_smooth(se = FALSE),
          'lm' = geom_smooth(se = FALSE, method = 'lm', formula = formula))
      }
    }
    if('ab_line' %in% input$plot_checkboxs) out <- out + geom_abline(slope = input$ab_slope, intercept = input$ab_intercept, size = 1, col = 'darkgrey', linetype = 2)
    out
    })
  output$plot_lm_res <- renderPrint({
    data_plt[,ratio := get(input$y_numerator)]
    if(input$y_denumerator != '1') data_plt[,ratio := ratio/get(input$y_denumerator)]
    
    x <- input$x
    if(input$polynom > 1){ for(i in 2:input$polynom) x <- paste0(x, ' + I(', input$x,'^', i, ')')}

    formula <- paste0('ratio ~ ', ifelse(input$intercept, '', '0 + '), x)
    formula <- as.formula(formula)
    out <- lm(formula, data_plt)

    formula_ <- paste0(input$y_numerator, ' / ', input$y_denumerator, ' ~ ', ifelse(input$intercept, '', '0 + '), x)
    out$call <- as.call(str2lang(formula_))
    out <- summary(out)
    out
  })
  output$freestyle_lm_res <- renderPrint({
    formula <- as.formula(input$lm_formula)
    out <- lm(formula, data_plt)
    out$call <- as.call(formula)
    out <- summary(out)
    out
    })
  output$data <- DT::renderDataTable(data)
}