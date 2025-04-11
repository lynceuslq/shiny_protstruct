#.libPaths('/mnt/mambaforge/envs/r_envs_shiny/lib/R/library')
.libPaths("/usr/local/lib/R/site-library")

library(shiny)
library(shinybusy)
library(bslib)
library(shinyBS)
library(shinycssloaders)
library(shinyjs)
library(shinyWidgets)
library(ggpubr)
library(survminer)
library(tidyverse)
library(rlang)
library(plyr)
library(tippy)
library(plotly)
library(kableExtra)
library(reactable)
library(NGLVieweR)
library(DiagrammeR)
library(ggmsa)
library(ggtree)
library(Biostrings)

# Define UI for application that draws a histogram
page_navbar(
  theme = bs_theme(version = 5, bootswatch = "lumen"),
  id="main",
  
  tags$style(
    ".warning-text {
      color: orange;
    }
    #element {
      color: orange;
    }
    "
  ),
  
  # Application title
  title="Visual Analysis of Protein Structure with Shiny",
  fillable="Dashboard",
  fillable_mobile=T,
  collapsible = T,
  nav_panel("Home",
            layout_column_wrap(
              width = 1,
              fill = FALSE,
              heights_equal = c("row"),
              h1("Welcome!"),
              h4(class="warning-text",icon("cov",class="fa-solid fa-circle-exclamation",lib = "font-awesome",style="font-size: 24px"),"Warning: the page is still under development."),
             DiagrammeROutput("workflow", height = "300px"),
              h4("Please click on the button and choose you analysis options to start your analysis."),
              
             uiOutput("start_project"),
             uiOutput("start_analysis")
             
            )
  ),
  nav_spacer(),
  
  nav_panel("Compare Structures",
            
            card(card_header("Compare uploaded structures"),
                 layout_sidebar(
                   fillable = TRUE),
                 full_screen = F)
            
  ),
  sidebar = sidebar(position = "right",
                    width = 300,
                    open = "closed",
                    h5("Need to report a problem?"),
                    h6("Developer Contacts:"),
                    h6(icon("cov",class="fa-regular fa-envelope",lib = "font-awesome",style="font-size: 24px"),
                       "qianli@zhejianglab.org"),
                    h6(icon("cov",class="fa-brands fa-github",lib = "font-awesome",style="font-size: 24px"),
                       "https://github.com/lynceuslq")
  )
  
)
