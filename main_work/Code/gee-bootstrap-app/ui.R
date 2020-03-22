library(shiny)
library(shinydashboard)

##### Header #####
header <- dashboardHeader()


##### Sidebar #####
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
    menuItem("Widgets", tabName = "widgets", icon = icon("th")),
    menuItem("Raw Data", tabName = "data", icon = icon("th"))
  ),
  selectInput('x', 'X Axis', sidebar_options[c(
    'P',
    'N',
    'NH',
    'ND',
    'P/N',
    'Length of Times Serie',
    'Percentage',
    'Actual SD'
    )]),
  selectInput('y_numerator', 'Y Axis Numerator', sidebar_options[c(
    'Actual SD',
    'GEE Estimated SD',
    'GEE(old) Estimated SD',
    'MLE Estimated SD'
  )]),
  selectInput('y_denumerator', 'Y Axis Denumerator', sidebar_options[c(
    '1',
    'Actual SD',
    'GEE Estimated SD',
    'GEE(old) Estimated SD',
    'MLE Estimated SD'
  )]),
  selectInput('color', 'Color', sidebar_options[c(
    'None',
    'P',
    'N',
    'NH',
    'ND',
    'P/N',
    'Length of Times Serie',
    'Percentage'
  )]),
  selectInput('aggFun', 'Trend', list(
    'None' = 'none', 'Mean' = 'mean', 'Median' = 'median', 'Smooth' = 'smooth', 'LM' = 'lm'
    )),
  radioButtons('vline', 'Vertical Line', list('None' = 0, '0.05' = 0.05, '1' = 1)),
  checkboxGroupInput(
    'plot_checkboxs', 'Plot Options',
    choices = list('Render Boxplot' = 'boxplot',
                   '0-1 Line' = 'abline',
                   'Jitter' = 'jitter'))
  )


##### Body #####
body <- dashboardBody(tabItems(
  tabItem(tabName = "dashboard", h2(
    plotOutput('graph'),
    box(
      verbatimTextOutput('plot_lm_res'),
      box(
        checkboxInput('intercept', 'Intercept', FALSE),
        numericInput('polynom', 'Polynom Degree', 1, 1, Inf, 1)
      ),
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
