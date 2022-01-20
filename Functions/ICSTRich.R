ICSTRich <- function(x, y){
  
  #X is a node Y is the total list of nodes. ICSRich is a function that returns the serotype richness 
  #Interclade serotype richness is calculated by determining how many serotypes are in all subclades
  
  #ICSRich is called and returns a character vector of serotypes
  
  # if( dim(TempData1)[1]==2 ){ 
  #Start at x
  TempData1 <- y[y$parent==x & y$node!=x , ]
  # BL <- c()
  if( dim(TempData1)[1]==2 ){
    
    if( all(is.na(TempData1$label) ) ){  #If double clade:
      BL <- c( ICSTRich(TempData1$node[1],y) , ICSTRich(TempData1$node[2],y))
      
      
    } else if( all(!is.na(TempData1$label)) ){  #If double tip
      
      BL <- c( as.character(TempData1$ST) )
      #
      
      
    } else { #Clade + tip
      
      BL <- c( ICSTRich( TempData1$node[ c(which(is.na(TempData1$label))) ],y) , as.character(TempData1$ST[ !is.na(TempData1$ST)   ]))
      
      
    }
    
    
    BL
    # } else {
    #   
    #   0
    # }
  }
  
}