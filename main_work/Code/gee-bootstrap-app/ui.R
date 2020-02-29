library(shiny)
library(shinydashboard)

##### Header #####
header <- dashboardHeader()


##### Sidebar #####
sidebar <- dashboardSidebar(
  selectInput('x', 'X Axis', sidebar_options[c('P', 'N', 'Percentage')]),
  selectInput('y', 'Y Axis', sidebar_options[c('Actual SD', 'GEE Estimated SD', 'Ratio')]),
  selectInput('color', 'Color', sidebar_options[c('None', 'P', 'N', 'Percentage')]),
  selectInput('shape', 'Shape', sidebar_options[c('None', 'P', 'N', 'Percentage')])
)


##### Body #####
body <- dashboardBody(
  plotOutput('plt'),
  dataTableOutput('data')
)


##### UI #####
dashboardPage(
  header,
  sidebar,
  body
)
