ICdist <- function(x, y){
  
  
  #X is a node Y is the total list of nodes. Branchlike is a function that returns the "interclade" distance 
  #Interclade distance is calculated by finding the mean length between tips
  
  
  #Start at x
  TempData1 <- y[y$parent==x , ]
  if( dim(TempData1)[1]==2 ){ 
    
    if( all(is.na(TempData1$label) ) ){  #If double clade:
      BL <- ( ICdist(TempData1$node[1],y)/2 + ICdist(TempData1$node[2],y)/2 + sum(TempData1$branch.length) )
      
      
    } else if( all(!is.na(TempData1$label)) ){  #If double tip
      
      BL <- sum(TempData1$branch.length)
      #
      
      
    } else { #Clade + tip
      
      BL <- ( ICdist( TempData1$node[ c(which(is.na(TempData1$label))) ],y)/2  + sum(TempData1$branch.length) ) 
      
      
    }
    
    
    BL
  } else {
    #Return 0 at root node, since this one is weird
    0
  }
  
}