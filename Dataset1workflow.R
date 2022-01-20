library(ggplot2)
library(ggtree)
library(treeio)
library(tidytree)
library(readr)
library(viridis)
library(ggnewscale)

Usedir <- dirname(rstudioapi::getSourceEditorContext()$path)

#Examples using dataset 1: 

Metadata <- readRDS(paste0(Usedir, "/Metadata/Dataset1meta.Rdata"))
Treefile <- paste0(Usedir, "/Metadata/Dataset1parsnp.tree")


Mytree <- read.newick(Treefile)

Metadata$Seqname <- gsub(".txt","",Metadata$Seqname)
Metadata$Seqname <- gsub("GCF","GCA",Metadata$Seqname)



x <- as_tibble(Mytree)
x$label <- gsub(".ref","",x$label)
d <- tibble(label = Metadata$Seqname, Serotype = Metadata$Serotype, ST = as.factor(Metadata$ST), 
)
y <- full_join(x, d, by = 'label')
# Mytree2 <- as.treedata(y)




# Tester <- as_tibble(Mytree2)


y$Serotype[which((y$Serotype)=="NT") ] <- NA

source(paste0(Usedir,"/Functions/ICdist.R"))
source(paste0(Usedir,"/Functions/ICSRich.R"))
source(paste0(Usedir,"/Functions/ICSTRich.R"))

#Define variables
y$ICD <- 0
y$ICSR <- 0
y$n <- 0
y$SD <- NA
y$SND <- NA

#Calculate 
for(i in 1:dim(y)[1]){
  
  if( is.na(y[i,]$label) ){
    
    
    
    y$ICD[i] <- ICdist(y$node[i],y)
    y$ICSR[i] <- length(levels(as.factor(ICSRich(y$node[i],y))))
    y$n[i] <- length(ICSRich(y$node[i],y))
    
    if( y$ICSR[i] >=2 ){
      y$SD[i] <- (y$ICSR[i]/y$ICD[i])
      y$SND[i] <- ((y$ICSR[i]/y$n[i])/y$ICD[i])
    }
    
    
  }
  
}


rm(i);
#All NaN to NA
y$ICSR[y$ICSR==0] <- NA
y$ICSR[y$ICSR==1] <- NA
y$n[y$n==0] <- NA
y$SD[which(is.nan(y$SD))] <- NA
y$SND[which(is.nan(y$SND))] <- NA
y$ICD[y$ICD==0] <- NA


#Reorder according to SND, descending 
Temp1 <- y[order(y$SND,decreasing = TRUE),]

Results <- data.frame()
for(i in 1:length(Temp1$SND[!is.na(Temp1$SND)])){
  
  Results <- rbind(Results, cbind(Temp1$node[i],
                                  paste0(levels(as.factor(ICSTRich(Temp1$node[i],Temp1))),collapse="."),
                                  paste0(levels(as.factor(ICSRich(Temp1$node[i],Temp1))),collapse="."), 
                                  Temp1$ICSR[i], 
                                  Temp1$n[i], 
                                  Temp1$SND[i]  ))
  
  
  
  
}
colnames(Results) <- c("Node","ST","Serotypes","ICSR","n","SND")


Mytree2 <- as.treedata(y)


Myserotype <- data.frame(as.factor(Metadata$Serotype))
colnames(Myserotype) <- "Serotype"
rownames(Myserotype) <- Metadata$Seqname

#Build circular tree, annotated with sequence type
p1 <- ggtree(Mytree2, layout="circular",branch.length = "none", size = .25) +
  geom_tiplab(size = 2.0, aes(x = 31-4, label=ST),check_overlap = TRUE) 

#Annotate with serotype
p1 <- gheatmap(p1, Myserotype, offset = 0, width = 0.1, colnames = FALSE, color = "Black") 

#
p1 <- p1 + ggnewscale::new_scale_fill() +
  geom_nodepoint(size = 5, shape = 21, color = "black", stroke = 0.3, aes(subset = (!is.na(SND))&SND>10,alpha = 0.6, fill = (SND) )) +  # geom_nodepoint(size = 0.7, stroke = 0.01, shape = 21, color = "black", fill = NA, aes(subset = (!is.na(SND))&SND>10 )) +
  scale_fill_viridis(alpha=0.8,option = "D", direction = -1,name="SND") +
  scale_alpha(guide="none")

p1
