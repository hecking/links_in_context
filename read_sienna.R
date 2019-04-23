require(igraph)
require(dplyr)

readGraphSienna <- function(file, mode="undirected") {
  
  read.table(file) %>% as.matrix %>% graph_from_adjacency_matrix(mode = mode) %>% 
    set_vertex_attr("name", value=1:vcount(.))
}

toDummyVars <- function(data) {
  
  lapply(1:ncol(data), function(col) {
    
    sapply(unique(data[,col]), function(val) {
      
      ifelse(data[,col] == val, 1, 0)
    })
  }) %>%
    do.call(cbind, .)
}

setDimNamesLazega <- function(atts) {
  
  rownames(atts) <- 1:nrow(atts) %>% as.character
  colnames(atts) <- c(
    "status=partner", "status=associate", "gender=male", "gender=female", 
    "office=Boston", "office=Hartford", "office=Providence",
    "years_in_firm=0-2", "years_in_firm=3-5", "years_in_firm=6-9", "years_in_firm=10-14", 
    "years_in_firm=15-20", "years_in_firm=21-25", "years_in_firm=26-30", "years_in_firm=31-35",
    "age=25-30", "age=31-35", "age=36-40", "age=41-45", "age=45-50", "age=51-55", "age=56-60", 
    "age=61-65", "age=66-71", "practice=litigation", "practice=corporate",
    "law_school=Harvard-Yale", "law_school=Ucon", "law_school=Other"
  )
  atts
}

readAttributesSienna <- function(file) {
  
  atts <- read.table(file) %>% select(-1) %>% 
    mutate(V5 = as.numeric(cut(V5, c(0,3,6,10,15,21,26,31,36)))) %>%
    mutate(V6 = as.numeric(cut(V6, c(25,31,36,41,46,51,56,61,66,71)))) %>% 
    toDummyVars %>%
    setDimNamesLazega
}