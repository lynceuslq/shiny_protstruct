#.libPaths('/mnt/mambaforge/envs/r_envs_shiny/lib/R/library')
.libPaths("/usr/local/lib/R/site-library")

library(shiny)
library(shinybusy)
library(bslib)
library(rlang)
library(plotly)
library(plyr)
library(dplyr)
library(foreach)
library(NGLVieweR)
library(ggmsa)
library(ggtree)
#library(readxl)
#library(xlsx)
#library(BiocManager)
#options(repos = BiocManager::repositories())
#ibrary(genbankr)
library(htmltools)
library(DT)
library(Biostrings)
library(foreach)
#library(DiagrammeR)
#library(DiagrammeRsvg)
#library(rsvg)
source("global.R")

function(input, output, session) {
  rv <- reactiveValues(query=NULL,
                       analysisdir=paste(Sys.Date(),format(Sys.time(), "%H-%M-%S"),sep="-"),
                       foldseek="/mnt/tools/foldseek/bin/foldseek",
                       foldmason="/mnt/tools/foldmason/bin/foldmason",
                       iqtree2="/mnt/tools/iqtree-2.3.6-Linux-intel/bin/iqtree2",
                       pdb_db="/mnt/pdb_db/pdb",
                       pdb2fasta="/mnt/tools/pdb2fasta",
                       interproscan="/mnt/tools/my_interproscan/interproscan-5.72-103.0/interproscan.sh",
                       outfile=NULL,
                       filter_fident=0,
                       filter_alntmscore=0,
                       filter_lddt=0,
                       num_target=0,
                       num_clu=0,
                       cluster_res=NULL,
                       tree_res=NULL,
                       iqtreedir=NULL,
                       targetdir=NULL)
  
  output$workflow <- renderDiagrammeR({
    
    return(diagramme)
  })
  
  observeEvent(input$goquery, {
    req(input$queryfiles)
   # rv$tmpdir <- paste("/mnt/shinyprocess/struct_tree",rv$analysisdir,sep = "/")
    rv$tmpdir <- paste("/home/vdb/shinyapps/shinyprocess/struct_tree",rv$analysisdir,sep = "/")
    files <- input$queryfiles
    system(paste("mkdir",rv$tmpdir))
    rv$refdir <- paste(rv$tmpdir,"Ref_structs",sep = "/")
    system(paste("mkdir",rv$refdir))
    print(paste0("storing analysis process in ", rv$tmpdir))
    for(nr in 1:length(files[, 1])){
      
      ext <- tools::file_ext(files[[nr, 'datapath']])
      print(ext)
     # shiny::validate(need(ext %in% c("pdb"), "Please upload files in .pdb format"))
      
      if(ext %in% c("pdb")) {
        targetfile <- paste0(rv$refdir,"/", files[[nr, 'name']])
        system(paste0("cp ",files[[nr, 'datapath']]," ", targetfile))
        rv$query <- c(rv$query,list(files[[nr, 'name']]))
      }
    }
    
    print(rv$query)
   # removeModal()
  })
  
  output$structureNGL_ref <- renderNGLVieweR({
   # req(rv$vis_ref)
   # req(rv$refdir)
    print(input$vis_ref)
    
    pdbpath <- paste0(rv$refdir,"/", input$vis_ref)
    print("to viz")
    print(pdbpath)
    
    ngl <-  NGLVieweR(pdbpath) %>%
      stageParameters(backgroundColor = "white") %>%
      setQuality("high") %>%
      setSpin(FALSE) %>%
      addRepresentation("cartoon",
                        param = list(
                          name = "cartoon",
                          colorScheme = input$colorScheme,
                          assembly=input$assembly
                        ) 
      ) %>%
      addRepresentation("label",
                        param = list(
                          name = "label", sele = "200:A.O",
                          showBackground = TRUE,
                          backgroundColor = "black",
                          backgroundMargin = 2,
                          backgroundOpacity = 0.5,
                          showBorder = TRUE,
                          colorValue = "white"
                        )
      )
    
    if (input$addsticks) {
      ngl <- ngl %>%
        addRepresentation("ball+stick",
                          param = list(
                            name ="ligand",
                            sele = "ligand",
                            colorValue = "red",
                            colorScheme = "element",
                            sele = "200"
                          )
        )
    }
    
    if(input$addsurface) {
      ngl <- ngl %>%
        addRepresentation("surface",
                          param = list(
                            name = "surface",
                            colorValue = "white",
                            opacity = 0.1
                          )
        ) 
    }
    
    if(!is.null(input$contact)) {
      newcontact <- set_contact
      newcontact[names(newcontact)==input$contact] <- T
      ngl <- ngl %>%
        addRepresentation("contact",
                          param = c(list(name = "contact"),newcontact)
        )
    }
    
    return(ngl)
  })
  
  output$start_analysis <- renderUI({
    if(is.null(rv$query)) {
      div()
    }else{
      div(layout_column_wrap(width=1/3,actionBttn("analysis_choice", "Choose analysis",color = "success"),
                             ),
          layout_column_wrap(width=1/3,actionBttn("clear_all", "Clear project",color = "warning"),
          ))
    }
  })
  
  observeEvent(input$analysis_choice, {
    showModal(choiceModal())
    
  })
  
  choiceModal <- function(failed = FALSE) {
    modalDialog(size="l",
                title="Strutural analysis",
                helpText("Please choose one of the analyses below:"
                ),
                layout_column_wrap(
                  width = 1/3,
                  fill = FALSE,
                  value_box(
                    "Search for similar strutures in PDBe",
                    #  uiOutput("total_tippers", container = h2),
                    #  textInput("depname","Plsease enter the name of analysis"),
                    actionBttn("ana_search","GO",color = "success"),
                    showcase = icon("cov",class="fa-solid fa-magnifying-glass",lib = "font-awesome",style="font-size: 24px")
                    #bsicons::bs_icon("bar-chart-steps")
                  ),
                  value_box(
                    "Structure clustering analysis",
                    actionBttn("ana_cluster","GO",color = "warning"),
                    showcase = icon("cov",class="fa-solid fa-circle-nodes",lib = "font-awesome",style="font-size: 24px")
                    #bsicons::bs_icon("bar-chart-steps")
                  ),
                  value_box(
                    "Tree inference based on structural alignments",
                    actionBttn("ana_tree","GO",color = "primary"),
                    showcase = icon("cov",class="fa-solid fa-sitemap",lib = "font-awesome",style="font-size: 24px")
                    #bsicons::bs_icon("bar-chart-steps")
                  )),
                  footer = tagList(
                    modalButton("Cancel")
                  #  actionButton("gomcmc", "Confirm")
                  )
    )
  }
  
  checkModal  <- function(failed = FALSE) {
    modalDialog(size="m",
                title="Check analysis",
                h4("Whoops! It seem this analysis has been conducted.
                   Click continue to remove previous data and start a new one."),
                actionButton("clearana", "Continue"))
  }
  
  observeEvent(input$clearana, {
    if(rv$oncheckana == "ana_search") {
      
      removeTab(inputId = "main", target = paste0("Structural alignments"))
      
      system(paste0("rm ", rv$outfile))
      
      if(rv$num_target > 0){
        system(paste0("rm -r ",rv$targetdir))
        system(paste0("rm ", rv$numoutfile))
      }
      
      rv$filter_fident=0
      rv$filter_alntmscore=0
      rv$filter_lddt=0
      rv$num_target=0
      
      rv$outfile <-  NULL
      rv$numoutfile <-  NULL
      rv$targetdir <- NULL
    }
    
    if(rv$oncheckana == "ana_cluster") {
      removeTab(inputId = "main", target = paste0("Structural clusters"))
      
      system(paste0("rm -r", rv$cluster_res))
      
      rv$num_clu <- 0
      rv$cluster_res <- NULL
    }
    
    removeModal()
  })
  
  warningModal <- function(failed = FALSE) {
    modalDialog(size="m",
                title="Are you sure to clear all data",
                actionButton("goclear", "Confirm"))
  }
  
  
  observeEvent(input$clear_all, {
    showModal(warningModal())
  })
  
  observeEvent(input$goclear, {
    rv$query <- NULL
    system(paste0("rm -r ",rv$tmpdir))
    rv$tmpdir <- NULL
    removeTab(inputId = "main", target = paste0("Query protein strutures"))
    removeTab(inputId = "main", target = paste0("Structural alignments"))
    removeTab(inputId = "main", target = paste0("Structural clusters"))
    
    rv$filter_fident=0
    rv$filter_alntmscore=0
    rv$filter_lddt=0
    rv$num_target=0
    rv$num_clu =0
    
    rv$cluster_res=NULL
    rv$outfile <-  NULL
    rv$numoutfile <-  NULL
    rv$targetdir <- NULL
    
    removeModal()
  })
  
  output$start_project <-renderUI({
    if(is.null(rv$tmpdir)) {
      div(
        layout_column_wrap(width=1/3,
                           actionBttn("start", "New Project")
        )
      )
    }else{
      div()
    }
  })
  
  observeEvent(input$goquery, {
    if(!is.null(rv$query)) {
      insertTab(inputId = "main",
                navbarMenu(paste0("Query protein strutures"),
                           
                           tabPanel("Struture visualization", 
                                    card(card_header("NGLVieweR"),
                                         layout_sidebar(
                                           fillable = TRUE,
                                           sidebar = sidebar(open=F,
                                                             selectInput("vis_ref",
                                                                         "Please select structure to visualize",
                                                                         choices = as.character(rv$query)),
                                                             selectInput("assembly",
                                                                         "Please select assembly to visualize",
                                                                         choices = c("AU","BU1","UNICELL","SUPERCELL")),
                                                             selectizeInput("colorScheme",
                                                                            "Please select color schemes",
                                                                            choices = c("atomindex","bfactor","chainid","chainindex","chainname",
                                                                                        "densityfit","electrostatic","element","entityindex",
                                                                                        "entitytype","geoquality","hydrophobicity","modelindex",
                                                                                        "moleculertype","occupancy","partialcharge"),
                                                                            ),
                                                             checkboxInput("addsurface","Molecular surface"),
                                                             checkboxInput("addsticks","Show ligand"),
                                                             checkboxGroupInput("contact","Show contact",
                                                                                choices = names(set_contact))),
                                           
                                           NGLVieweROutput("structureNGL_ref")
                                         ),
                                         full_screen = T)
                                    
                           ))
      )
    }
    
    
    removeModal()
  })
  
  
  dataModal <- function(failed = FALSE) {
    modalDialog(size="l",title="Query strutures",
                helpText("Please upload your protein strutures for analysis (. file only):"
                ),
                fileInput('queryfiles', 
                          "Query proteins:",
                          multiple = T,
                          accept = c(".pdb")),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("goquery", "Confirm")
                )
    )
  }
  
  observeEvent(input$start, {
    showModal(dataModal())
  })
  
  
  #### structure search analysis
  observeEvent(input$ana_search, {
    rv$oncheckana <- "ana_search"
    if(is.null(rv$outfile)) {
      showModal(searchparsModal())
    }else{
      showModal(checkModal())
    }
  })
  
  
  searchparsModal <- function(failed = FALSE) {
    modalDialog(size="l",
                title="Choose paramters for  search with Foldseek",
                searchpar_div,
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("gosearch", "Confirm")
                )
    ) 
  }
  
  observeEvent(input$gosearch, {
    rv$refdir <- paste(rv$tmpdir,"Ref_structs",sep = "/")
    rv$outfile <-  paste(rv$tmpdir,"struct.aln.m8",sep = "/")
    rv$numoutfile <-  paste(rv$tmpdir,"num_searchout.list",sep = "/")
    tmpFolder <- paste(rv$tmpdir,"tmpFolder",sep = "/")
    show_modal_spinner() 
    system(paste(
      rv$foldseek,
      "easy-search",
      rv$refdir,
      rv$pdb_db,
      "-s",
      input$search_sens,
      "--min-seq-id",
      input$search_min_id,
      "--max-seqs",
      input$search_max_seq,
      " --format-output 'query,target,fident,alnlen,alntmscore,qtmscore,ttmscore,lddt,prob,qstart,qend,tstart,tend,evalue,bits'",
      rv$outfile,
      tmpFolder
      )
    )
    remove_modal_spinner() 
    print("easy-search completed")
    system(paste("cat",rv$outfile, " | wc -l > " ,rv$numoutfile ))
    rv$num_out <- read_lines(rv$numoutfile)
   # rv$num_out <- as.numeric(try(system(paste("cat",rv$outfile, " | wc -l"))))
    print(paste("number of hits",rv$num_out))
    showModal(searchcheckModal())
  })
  
  searchcheckModal <- function(failed = FALSE) {
    modalDialog(size="m",
                title=paste(rv$num_out, "hits were found!"),
                actionButton("gonext", "Continue"))
  }
  
  rv$targettab <- reactive({
    req(input$gonext)
    req(rv$num_out)
    if(!is.null(rv$query)) {
      if(rv$num_out > 0) {
        alntab <- read.table( rv$outfile,sep="\t")
        colnames(alntab) <- c("query","target","fident","alnlen","alntmscore","qtmscore",
                              "ttmscore","lddt","prob","qstart","qend","tstart","tend","evalue","bits")
       # print(alntab)
        alntab <- subset(alntab,fident>=rv$filter_fident & alntmscore>= rv$filter_alntmscore & lddt >= rv$filter_lddt)
      return(alntab)
        }}
  })
  
  rv$selectedtargets <- reactive({
    req(rv$targettab())
    tartab <- rv$targettab()
    if(length(tartab$target) >0) {
      targets <- unique(tartab$target)
    #  print(paste("targets",targets))
      uppertagets <- foreach::foreach(a=targets,.combine = "c") %do% toupper(strsplit(a,"[.]")[[1]][1])
      print(unique(uppertagets))
      
      return(unique(uppertagets))
    }else{
      return(NULL)
    }
    
  })
  
  output$aln_tab <- renderReactable({
    req(input$gonext)
      print("reactable")
      print(rv$outfile)
      if(!is.null(rv$query)) {
        if(rv$num_out > 0) {
        
        print(rv$targettab())
        reactable(rv$targettab(),
                  groupBy = c("query"),
                  columns = list(
          target = colDef(cell = function(value, index) {
            # Render as a link
            id=toupper(strsplit(value,"[.]")[[1]][1])
            url <- paste0("https://www.rcsb.org/structure/", id)
            htmltools::tags$a(href = url, target = "_blank", as.character(id))
          }),
          alntmscore = colDef(
            style = function(value) {
              if (value > 0.5) {
                color <- "#008000"
              } else if (value < 0.3) {
                color <- "#e00000"
              } else {
                color <- "#777"
              }
              list(color = color, fontWeight = "bold")
            }
          ),
          fident = colDef(
            style = function(value) {
              if (value > 0.5) {
                color <- "#008000"
              } else if (value < 0.3) {
                color <- "#e00000"
              } else {
                color <- "#777"
              }
              list(color = color, fontWeight = "bold")
            }
          ),
          lddt = colDef(
            style = function(value) {
              if (value > 0.5) {
                color <- "#008000"
              } else if (value < 0.3) {
                color <- "#e00000"
              } else {
                color <- "#777"
              }
              list(color = color, fontWeight = "bold")
            }
          )
          )
        )
        }}
      
  })
  
  observeEvent(input$gonext, {
    if(!is.null(rv$query)) {
      if(rv$num_out > 0) {
        insertTab(inputId = "main",
                  navbarMenu(paste0("Search Results"),
                             
                             tabPanel("Structure aligned to query proteins", 
                                      card(card_header("Hits"),
                                           layout_column_wrap(width=1/3,
                                                              actionBttn("filtertarget", "Filter Targets"),
                                                              actionBttn("clearfilter", "Clear Filtering"),
                                                              actionBttn("selecttargets", "Confirm Selection")),
                                           reactableOutput("aln_tab"),
                                           full_screen = T)
                                      
                             ))
        )
      }
      
    }
    
    
    removeModal()
  })
  
  observeEvent(input$filtertarget, {
    showModal(filtertargetModal())
  })
  
  filtertargetModal <- function(failed = FALSE) {
    modalDialog(size="m",
                title="Filter structural alignment hits",
                sliderInput("filter_fident","Choose targets with fident score higher than",
                            min=0,
                            max=1,
                            value=0.5,
                            step =0.01),
                sliderInput("filter_alntmscore","Choose targets with alntmscore score higher than",
                            min=0,
                            max=1,
                            value=0.5,
                            step =0.01),
                
                sliderInput("filter_lddt","Choose targets with lddt score higher than",
                            min=0,
                            max=1,
                            value=0.5,
                            step =0.01),
                
                actionButton("gotargetfilter", "Confirm"))
  }
  
  observeEvent(input$selecttargets, {
    showModal(selecttargetsModal())
  })
  selecttargetsModal <- function(failed = FALSE) {
    modalDialog(size="m",
                title = paste(length(rv$selectedtargets()),"targets were selected!"),
                h4(paste("Click on continue to download target protein structures")),
                actionButton("complete_aln", "Continue"))
  }
  
  observeEvent(input$complete_aln, {
    rv$targetdir <- paste(rv$tmpdir,"target_structs",sep = "/")
    system(paste("mkdir",rv$targetdir))
    
    targets <- rv$selectedtargets()
    if(!is.null(targets) ) {
      rv$num_target <- length(targets)
      show_modal_progress_line(text = "Start downloading")
      for(t in 1:length(targets)) {
        system(paste0("cd ",
                      rv$targetdir,
                      " ; ",
                      " wget https://files.rcsb.org/download/",
                      targets[t],
                      ".pdb ; echo -e ",
                      targets[t])
        )
        update_modal_progress(
          value = t / length(targets)
        )
      }
      remove_modal_progress()
    }
    
    removeModal()
  })
  
  observeEvent(input$gotargetfilter,{
      rv$filter_fident <- input$filter_fident
      rv$filter_alntmscore <- input$filter_alntmscore
      rv$filter_lddt <- input$filter_lddt
    
    
    print(paste("filtering at:",
                "fident:", rv$filter_fident,
                "alntmscore:", rv$filter_alntmscore,
                "lddt:", rv$filter_lddt
    ))
    
    removeModal()
    
  })
  
  observeEvent(input$clearfilter,{
        rv$filter_fident <- 0
        rv$filter_alntmscore <- 0
        rv$filter_lddt <- 0
      
    })
  
  ### structure clustering analysis

  observeEvent(input$ana_cluster, {
    rv$oncheckana <- "ana_cluster"
    if(is.null(rv$cluster_res)) {
      showModal(clusterparModal())
    }else{
      showModal(checkModal())
    }
  })
  
  clusterparModal <- function(failed = FALSE) {
    if(!is.null(rv$targetdir)) {
      modalDialog(size="m",
                  title = "Set parameter for srtructure clustering",
                  checkboxInput("include_target","Include selected targets for clustering",value=T),
                  cluster_pars_div,
                  actionButton("com_clu_pars", "Continue"))
    }else{
      modalDialog(size="m",
                  title = "Set parameter for srtructure clustering",
                  cluster_pars_div,
                  actionButton("com_clu_pars", "Continue"))
    }
    
  }
 
  observeEvent(input$com_clu_pars,{
    tmpFolder <- paste(rv$tmpdir,"clustertmpFolder",sep = "/")
    rv$cluster_res <- paste(rv$tmpdir,"cluster_res",sep = "/")
    system(paste("mkdir",
                 rv$cluster_res))
    
    if(!is.null(rv$targetdir)) {
      
      if((rv$num_target + length(rv$query)) < 3) {
        showModal(noclusteringModal())
      }else{
        show_modal_spinner()
        system(paste0("cp ",
                     rv$refdir,
                     "/*.pdb ",
                     rv$targetdir))
        system(paste(
          rv$foldseek,
          "easy-cluster",
          rv$targetdir,
          paste0(rv$cluster_res,"/res"),
          tmpFolder,
          "--cov-mode",
          input$cov_mode,
          "-c",
          input$cov_fract,
          "--cluster-mode",
          input$cluster_mode,
          "--min-seq-id",
          input$clu_min_seq_id,
          "--tmscore-threshold",
          input$clu_tmscore_threshold,
          "--lddt-threshold",
          input$clu_lddt_threshold
        ))
        remove_modal_spinner() 
        rv$num_clu <- as.numeric(try(system(paste("cat",
                                                  paste(rv$cluster_res, "res_cluster.tsv", sep = "/" ),
                                                  " | wc -l "
        ),intern = TRUE)))
        showModal(clustercheckModal())
        print("number of clusters")
        print(rv$num_clu)
      }
    }else{
      if(length(rv$query) < 3) {
        showModal(noclusteringModal())
      }else{
        show_modal_spinner() 
        system(paste(
          rv$foldseek,
          "easy-cluster",
          rv$refdir,
          paste0(rv$cluster_res,"/res"),
          tmpFolder,
          "--cov-mode",
          input$cov_mode,
          "-c",
          input$cov_fract,
          "--cluster-mode",
          input$cluster_mode,
          "--min-seq-id",
          input$clu_min_seq_id,
          "--tmscore-threshold",
          input$clu_tmscore_threshold,
          "--lddt-threshold",
          input$clu_lddt_threshold
          ))
        
        remove_modal_spinner() 
        
        rv$num_clu <- as.numeric(try(system(paste("cat",
                                                  paste(rv$cluster_res, "res_cluster.tsv", sep = "/" ),
                                                  " |cut -f1 | sort | uniq | wc -l "
                                                  ),intern = TRUE)))
        
        
        showModal(clustercheckModal())
        print("number of clusters")
        print(rv$num_clu)
      }
    }
    
    
  })
  
  clustercheckModal <- function(failed = FALSE) {
    modalDialog(size="m",
                title=paste(rv$num_clu, "struture clusters were found!"),
                actionButton("gonextclu", "Continue"))
  }
  
  rv$clustermatch <- reactive({
    req(input$gonextclu)
    req(rv$num_clu)
    if(!is.null(rv$query)) {
        clustermatch <- read.table(paste(rv$cluster_res, "res_cluster.tsv", sep = "/" ),sep = "\t")
        colnames(clustermatch) <- c("Representative","Proteins")
        return(clustermatch)
      }
  })
  
  output$cluster_tab <- renderReactable({
    req(input$gonextclu)
    print(rv$clustermatch())
    mytab <- rv$clustermatch()
    
    if(!is.null(rv$query)) {
      print(as.character(rv$query))
        reactable(mytab,
                  groupBy = c("Representative"),
                  columns = list(
                    Proteins = colDef(
                      style = function(value) {
                        if (paste(value,"pdb",sep = ".") %in% rv$query) {
                          color <- "#008000"
                        } else {
                          color <- "#777"
                        }
                        list(color = color)
                      }
                    )
                  ))
    }
  })
  
  observeEvent(input$gonextclu, {
    if(!is.null(rv$query)) {
        insertTab(inputId = "main",
                  navbarMenu(paste0("Structural clusters"),
                             
                             tabPanel("Structure clustering results", 
                                      card(card_header("Clusters"),
                                           reactableOutput("cluster_tab"),
                                           full_screen = T)
                                      
                             )))
    }
    removeModal()
  })
  
  noclusteringModal <- function(failed = FALSE) {
    modalDialog(size="m",
                title = "Cannot cluster less than three protein structures!")
  }

#### alignment and tree
  observeEvent(input$ana_tree, {
    rv$oncheckana <- "ana_tree"
    showModal(treeparsModal())
  })
  
  treeparsModal <- function(failed = FALSE) {
    if(!is.null(rv$targetdir)) {
      modalDialog(size="m",
                  title = "Set parameter for alignments based on structures",
                  checkboxInput("tree_include_target","Include selected targets for alignments",value=T),
                  treepars_div,
                  actionButton("com_tree_pars", "Continue"))
    }else{
      modalDialog(size="m",
                  title = "Set parameter for alignments based on structures",
                  treepars_div,
                  actionButton("com_tree_pars", "Continue"))
    } 
  }
  
  observeEvent(input$com_tree_pars,{
    tmpFolder <- paste(rv$tmpdir,"treetmpFolder",sep = "/")
    rv$tree_res <- paste(rv$tmpdir,"tree_res",sep = "/")
    system(paste("mkdir", rv$tree_res))
    
    if(input$recompute_score) {
      aln_recompute_score <- 1
    }else{
      aln_recompute_score <- 0
    }
    
    if(input$regressive) {
      aln_regressive <- 1
    }else{
      aln_regressive <- 0
    }
    
    if(!is.null(rv$targetdir)) {
      
      if((rv$num_target + length(rv$query)) < 3) {
        showModal(nomaketreeModal())
      }else{
        show_modal_spinner()
        system(paste("mkdir",
                     paste(rv$tree_res,"structs",sep = "/")))
        system(paste0("cp ",
                      rv$refdir,
                      "/*.pdb ",
                      paste(rv$tree_res,"structs",sep = "/")))
        system(paste0("cp ",
                      rv$targetdir,
                      "/*.pdb ",
                      paste(rv$tree_res,"structs",sep = "/")))
        
        system(paste(
          rv$foldmason,
          "easy-msa",
          paste(rv$tree_res,"structs",sep = "/"),
          paste0(rv$tree_res,"/result.fasta"),
          tmpFolder,
          "--match-ratio",
          input$match_ratio,
          "--mask-bfactor-threshold",
          input$mask_bfactor_threshold,
          "--refine-iters",
          input$refine_iters,
          "--recompute-scores",
          aln_recompute_score,
          "--regressive",
          aln_regressive
        ))
        remove_modal_spinner() 
        showModal(treecheckModal())
        print("aln completed")
        rv$protein_sequences <- paste0(rv$tree_res,"/result.fasta_aa.fa")
        rv$aln_3di <- paste0(rv$tree_res,"/result.fasta_3di.fa")
        rv$AAMultipleAlignment <- readAAMultipleAlignment(rv$protein_sequences)
        rv$AAMultipleAlignment_3di <- readAAMultipleAlignment(rv$aln_3di)
        
        rv$aln_width <- as.data.frame(rv$AAMultipleAlignment@unmasked@ranges)[1,"width"]
        rv$aln_width_3di <- as.data.frame(rv$AAMultipleAlignment_3di@unmasked@ranges)[1,"width"]
        
      }
    }else{
      if(length(rv$query) < 3) {
        showModal(nomaketreeModal())
      }else{
        show_modal_spinner() 
        system(paste(
          rv$foldmason,
          "easy-msa",
          rv$refdir,
          paste0(rv$tree_res,"/result.fasta"),
          tmpFolder,
          "--match-ratio",
          input$match_ratio,
          "--mask-bfactor-threshold",
          input$mask_bfactor_threshold,
          "--refine-iters",
          input$refine_iters,
          "--recompute-scores",
          aln_recompute_score,
          "--regressive",
          aln_regressive
        ))
        remove_modal_spinner() 
        showModal(treecheckModal())
        rv$protein_sequences <- paste0(rv$tree_res,"/result.fasta_aa.fa")
        rv$aln_3di <- paste0(rv$tree_res,"/result.fasta_3di.fa")
        rv$AAMultipleAlignment <- readAAMultipleAlignment(rv$protein_sequences)
        rv$AAMultipleAlignment_3di <- readAAMultipleAlignment(rv$aln_3di)
        
        rv$aln_width <- as.data.frame(rv$AAMultipleAlignment@unmasked@ranges)[1,"width"]
        rv$aln_width_3di <- as.data.frame(rv$AAMultipleAlignment_3di@unmasked@ranges)[1,"width"]
        print("aln completed")
      }
    }
  })
  
  treecheckModal <- function(failed = FALSE) {
    modalDialog(size="m",
                title=paste("Alignment completed"),
                actionButton("gonexttree", "Continue"))
  }
  
  observeEvent(input$gonexttree, {
    if(!is.null(rv$query)) {
      if(!is.null(rv$tree_res)) {
      insertTab(inputId = "main",
                navbarMenu(paste0("Alignments"),
                           
                           tabPanel("Alignments based on structures", 
                                    layout_column_wrap(width=1/2,
                                                       card(card_header("Guide tree and alignments"),
                                                            layout_sidebar(
                                                              fillable = TRUE,
                                                              sidebar = sidebar(open=F,
                                                                                selectInput("showalnpl","Show alignments by",c("Amino Acid","3Di")),
                                                                                checkboxInput("gtree_cir","Show in circular layout")),
                                                              
                                                              plotOutput("gtree_plot"),
                                                              actionButton("refinetree", "Refine tree inference")
                                                            ),
                                                            full_screen = T),
                                                       card(card_header("Alignments Details"),
                                                            layout_sidebar(
                                                              fillable = TRUE,
                                                              sidebar = sidebar(open=F,
                                                                                selectInput("showaln","Show alignments by",c("Amino Acid","3Di")),
                                                                                sliderInput("aln_loci", label = "Showing loci", min = 1, 
                                                                                            max = rv$aln_width, value = c(1, min(50,rv$aln_width))),
                                                                                checkboxInput("aln_logo","Show sequence logo")),
                                                              
                                                              plotOutput("aln_plot")
                                                            ),
                                                            full_screen = T)
                                                       
                                    ),
                                    uiOutput("refinetree_div")
                           )
                           ))
    }}
    removeModal()
  })
  
  output$gtree_plot <- renderPlot({
    if(!is.null(rv$tree_res)) {
      if(input$showalnpl == "3Di") {
        sequences <- rv$aln_3di
      }else{
        sequences <- rv$protein_sequences
      }
      
      treefile <- paste0(rv$tree_res,"/result.fasta.nw")
      if(input$gtree_cir) {
        p <- ggtree(read.tree(text=paste0(read_lines(treefile),";")), layout='circular') + geom_tiplab(size=3)
        plot <- msaplot(p, sequences, window=c(120, 200))
      }else{
        p <- ggtree(read.tree(text=paste0(read_lines(treefile),";"))) + geom_tiplab(size=3)
        plot <- msaplot(p, sequences, offset=3, width=2)
      }
      return(plot)
    }
  })
  
  output$aln_plot <- renderPlot({
    if(!is.null(rv$tree_res)) {
      if(input$showaln == "3Di") {
        sequences <- rv$aln_3di
      }else{
        sequences <- rv$protein_sequences
      }
      
      if(input$aln_logo) {
        plot <-  ggmsa(sequences, start = input$aln_loci[1], end = input$aln_loci[2], char_width = 0.5, seq_name = T) + 
          geom_seqlogo()
      }else{
        plot <-  ggmsa(sequences, start = input$aln_loci[1], end = input$aln_loci[2], char_width = 0.5, seq_name = T) 
      }
      
      return(plot)
    }
  })
  
  observeEvent(input$refinetree,{
    showModal(refinetreeparsModal())
  })
  
  output$refinetree_div <- renderUI({
    if(!is.null(rv$iqtreedir)){
      
      if(rv$num_clu >0){
        card(card_header("Maximum likelihood tree"),
             layout_sidebar(
               fillable = TRUE,
               sidebar = sidebar(
                 selectInput("refinetreelayout","Select layout",
                             choices = c('rectangular', 'dendrogram', 'slanted',
                                         'ellipse', 'roundrect', 'fan', 'circular', 'inward_circular',
                                         'radial', 'equal_angle', 'daylight', 'ape')),
                 checkboxInput("includeclures","Color tree tips by clustering results")
               ),
               plotOutput("refinetree_plot")
             ),
             full_screen = T)
      }else{
        card(card_header("Maximum likelihood tree"),
             layout_sidebar(
               fillable = TRUE,
               sidebar = sidebar(
                 selectInput("refinetreelayout","Select layout",
                             choices = c('rectangular', 'dendrogram', 'slanted',
                                         'ellipse', 'roundrect', 'fan', 'circular', 'inward_circular',
                                         'radial', 'equal_angle', 'daylight', 'ape'))
               ),
               plotOutput("refinetree_plot")
             ),
             full_screen = T)
      }
    }else{
      div()
    }
  })
  
  output$refinetree_plot <- renderPlot({
    if(!is.null(rv$iqtreedir)) {
      treefile <- paste0(paste0(rv$iqtreedir,"/tree"),".treefile")
      if(rv$num_clu >0) {
        if(input$includeclures){
          mytab <- rv$clustermatch()
          
          grouplist <- foreach::foreach(a=unique(mytab$Representative)) %do% {
            mytab$Proteins[mytab$Representative==a]
          }
          names(grouplist) <- paste("cluster",1:length(unique(mytab$Representative)))
          
          print(grouplist)
          t <- ggtree(read.tree(text=read_lines(treefile)),
                      layout = input$refinetreelayout) +
            geom_tiplab(size=3)
          p <- groupOTU(t, grouplist, 'Group') + aes(color=Group) + 
            theme(legend.position = "none")
          
        }else{
          p <- ggtree(read.tree(text=read_lines(treefile)),
                      layout = input$refinetreelayout) +
            geom_tiplab(size=3)
        }
      }else{
        p <- ggtree(read.tree(text=read_lines(treefile)),
                    layout = input$refinetreelayout) +
          geom_tiplab(size=3)
      }
      return(p)
    }})
  
  refinetreeparsModal <- function(failed = FALSE) {
    
    modalDialog(size="m",
                  title = "Start to reconstruct maximum likelihood tree with iqtree2",
                  actionButton("refine_tree_pars", "Continue"))
  }
  
  observeEvent(input$refine_tree_pars,{
    if(!is.null(rv$tree_res)) {
      rv$iqtreedir <- paste0(rv$tree_res,"/iqtree")
      show_modal_spinner() 
      system(paste("mkdir",rv$iqtreedir))
    system(paste(
      rv$iqtree2,
      "-s",
      paste0(rv$tree_res,"/result.fasta_aa.fa"),
      "--seqtype",
      "AA",
      "-T AUTO",
      "--prefix",
      paste0(rv$iqtreedir,"/tree")
    ))
    remove_modal_spinner()
    print("iqtree completed")
    }
    
  })
  
  nomaketreeModal <- function(failed = FALSE) {
    modalDialog(size="m",
                title = "Cannot generate MSA with less than three proteins!")
  }
}