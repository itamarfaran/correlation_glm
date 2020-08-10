library(shiny)
library(shinydashboard)

##### Header #####
header <- dashboardHeader(title = 'Simulations Results')


##### Sidebar #####
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
    # menuItem("Widgets", tabName = "widgets", icon = icon("th")),
    menuItem("Raw Data", tabName = "data", icon = icon("th"))
  ),
  selectInput('x', 'X Axis', sidebar_options[c(
    'P',
    'M',
    'N',
    'NH',
    'ND',
    'AR Coefficient',
    'P/N',
    'M/N',
    'Length of Times Serie',
    'Percentage of Sick Samples',
    'None Null Alphas',
    'Minimum Range of Alpha',
    'Real Alpha',
    'Alpha Estimate',
    'Alpha Estimate SD'
    )]),
  selectInput('y_numerator', 'Y Axis Numerator', sidebar_options[c(
    'Alpha Estimate',
    'Alpha Estimate SD',
    'GEE Estimated SD',
    'SD of GEE Estimated SD'
    )]),
  selectInput('y_denumerator', 'Y Axis Denumerator', sidebar_options[c(
    '1',
    'Alpha Estimate',
    'Real Alpha',
    'Alpha Estimate SD',
    'GEE Estimated SD'
  )]),
  selectInput('color', 'Color', sidebar_options[c(
    'None',
    'P',
    'M',
    'N',
    'NH',
    'ND',
    'AR Coefficient',
    'P/N',
    'M/N',
    'Length of Times Serie',
    'Percentage of Sick Samples',
    'None Null Alphas',
    'Minimum Range of Alpha',
    'Real Alpha',
    'Alpha Estimate'
    )]),
  checkboxGroupInput(
    'plot_checkboxs', 'Plot Options',
    choices = list(
      'Color as Factor' = 'color_as_factor',
      'Render Boxplot' = 'boxplot',
      'Jitter' = 'jitter',
      'Add AB Line' = 'ab_line',
      'Drop Line at 0' = 'no_zero_line')
    ),
  splitLayout(
    cellWidths = c('50%', '50%'),
    numericInput('ab_intercept', 'AB Intercept', value = 0),
    numericInput('ab_slope', 'AB Slope', value = 0)
    ),
  selectInput('aggFun', 'Trend', list(
    'None' = 'none', 'Mean' = 'mean', 'Median' = 'median', 'Smooth' = 'smooth', 'LM' = 'lm'
  ))#,
  # selectInput('x_scale', 'X Axis Scale', names(scales_x)),
  # selectInput('y_scale', 'Y Axis Scale', names(scales_y))
)


##### Body #####
body <- dashboardBody(tabItems(
  tabItem(tabName = "dashboard", h2(
    plotOutput('graph'),
    box(
      splitLayout(
        checkboxInput('intercept', 'Intercept', FALSE),
        numericInput('polynom', 'Polynom Degree', 1, 1, Inf, 1)
      ),
      verbatimTextOutput('plot_lm_res'),
      title = 'LM on Plot',
      collapsible = TRUE
    ),
    box(
      textInput('lm_formula', 'Free Style LM Formula', 'gee_sd_mean/alpha_est_emp_sd ~ p + n'),
      verbatimTextOutput('freestyle_lm_res'),
      collapsible = TRUE
    )
  )),
  tabItem(tabName = "data", h2(
    dataTableOutput('data')
  ))
))

##### UI #####
dashboardPage(
  header,
  sidebar,
  body
)
