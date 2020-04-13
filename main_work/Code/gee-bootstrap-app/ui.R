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
    'P/N',
    'M/N',
    'AR Coefficient',
    'Length of Times Serie',
    'NH',
    'ND',
    'Percentage of Sick Samples',
    'Actual SD'
    )]),
  selectInput('y_numerator', 'Y Axis Numerator', sidebar_options[c(
    'Actual SD',
    'GEE Estimated SD'
    )]),
  selectInput('y_denumerator', 'Y Axis Denumerator', sidebar_options[c(
    '1',
    'Actual SD',
    'GEE Estimated SD'
  )]),
  selectInput('color', 'Color', sidebar_options[c(
    'None',
    'P',
    'N',
    'NH',
    'ND',
    'P/N',
    'AR Coefficient',
    'Length of Times Serie',
    'Percentage of Sick Samples'
  )]),
  checkboxGroupInput(
    'plot_checkboxs', 'Plot Options',
    choices = list(
      'Color as Factor' = 'color_as_factor',
      'Render Boxplot' = 'boxplot',
      'Jitter' = 'jitter',
      'Add AB Line' = 'ab_line')
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
      textInput('lm_formula', 'Free Style LM Formula', 'gee_sd/actual_sd ~ p + n'),
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
