set_contact <- list("hydrogenBond"=F,
                    "weakHydrogenBond"=F,
                    "waterHydrogenBond"=F,
                    "backboneHydrogenBond"=F,
                    "hydrophobic"=F,
                    "halogenBond"=F,
                    "ionicInteraction"=F,
                    "metalCoordination"=F,
                    "cationPi"=F,
                    "piStacking"=F,
                    "refineSaltBridges"=F)

cluster_pars_div <- div(
  
  selectInput("cov_mode","select covering mode",
              choices = c( "coverage of query and target"=0,
                           "coverage of target"=1,
                           "coverage of query"=2,
                           "target seq. length has to be at least x% of query length"=3,
                           "query seq. length has to be at least x% of target length"=4,
                           "short seq. needs to be at least x% of the other seq. length"=5),
              selected=1),
  sliderInput("cov_fract","Fraction of covered residues",
              min=0,max=1,
              value=0.6,
              step =0.01),
  
  selectInput("cluster_mode","Cluster mode",
              choices = c("greedy"=0,
                          "BLASTclust"=1,
                          "CDHIT"=2),
              selected=1),
  sliderInput("clu_min_seq_id","Alignment identity for clustering",
              min=0,max=1,
              value=0.5,
              step =0.01),
  sliderInput("clu_tmscore_threshold","Cluster alignments with a tmsore threshold",
              min=0,max=1,
              value=0.5,
              step =0.01),
  sliderInput("clu_lddt_threshold","Cluster alignments with a LDDT threshold",
              min=0,max=1,
              value=0.5,
              step =0.01)
)

searchpar_div <- div(
  sliderInput("search_sens", "Sensitivity",min=1,max=10,value=7.5,step = 0.5),
  sliderInput("search_min_id", "Minimum identity for hits",min=0.0,max=1.0,value=0.1),
  numericInput("search_max_seq", "Maximum number of hits per query",min=1,max=1000,value=10),
  
)

treepars_div <- div(
  sliderInput("match_ratio", "Match ratio",min=0,max=1,value=0.51,step = 0.01),
  sliderInput("mask_bfactor_threshold", "Mask bfactor threshold",min=0,max=100,value=0.00,step = 0.01),
  numericInput("refine_iters","refine iters",min=0,max=10000,value=0),
  checkboxInput("recompute_score","Recompute all-vs-all alignment scores every iteration"),
  checkboxInput("regressive","Align sequences root-to-leaf")
  
)




diagramme <- DiagrammeR::mermaid("
graph LR
  A{Home}-->V{New Project}
  V-->B[Query Proteins]
  B-->C{Select Analysis}
  B-->D(Strutural Visualization)
  C-->F(Structure Search)
  F-.->G[Select Targets]
  G-->C
  C-->J(Structure Clustering)
  J-.->K(Structure Alignment)
  C-->K(Structure Alignment)
  K-->L(Tree Reconstruction)
  
  J-->W{Clear Project}
  K-->W
  L-->W

  W-.->V
  
  
  style F fill:#ededed
  style G fill:#ededed
  style J fill:#ededed
  style K fill:#ededed
  style L fill:#ededed
")