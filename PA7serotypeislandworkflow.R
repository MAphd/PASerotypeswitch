
library(readr)

Usedir <- dirname(rstudioapi::getSourceEditorContext()$path)

Metadata1 <- readRDS(paste0(Usedir, "/Metadata/Dataset1meta.Rdata"))
Metadata2 <- readRDS(paste0(Usedir, "/Metadata/Dataset2meta.Rdata"))

source(paste0(Usedir,"/Functions/GIfinder.R"))




#Combine both datasets
Metadata1$Dataset <- "1"
Metadata2$Dataset <- "2"
Datafinal <- rbind(Metadata1, Metadata2)

#Only use O12 isolates subset
Datafinal <- subset(Datafinal, Datafinal$Serotype=="O12")




#GIfinder computes coordinates of O12 pa7-like serotype islands. If parameters Subject and Querydir are specified it will run blast alignments and 
#generate necessary data. 

dOutput <- paste0(Usedir, "/Alignments/")

#Minimum input
GIlist <- GIfinder(Outputdir = dOutput, Blasted = TRUE)

# dSubject <- paste0(Usedir,"/GCF_000017205.1_ASM1720v1_genomic.fna")

#Querydir is directory containing query O12 genomes, outputdir is where the resulting alignments will be saved. 
#Ensure that the outputdir is empty, since the script will attempt to read any .txt files in this directory. 

#GIlist <- GIfinder(Subject = paste0(Usedir,"/GCF_000017205.1_ASM1720v1_genomic.fna"), Querydir =, Outputdir =, Blasted = FALSE)


#Generate table of PA7 serotype island coordinates & pair with names
dMat <- data.frame()
for(i in 1:length(GIlist)){
  
  Tempvar <- as.numeric(GIlist[[i]][[1]][3:4])
  bleh <- which(  grepl(names(GIlist)[i], Datafinal$Seqname,fixed=TRUE)  ) 
  
  
  # plotobj <- plotobj + geom_segment(aes(x = Tempvar[1], xend=Tempvar[2],y=seq(1.2,1.2+(length(dInput)*0.1),0.1)[i], yend = seq(1.2,1.2+(length(dInput)*0.1),0.1)[i],size = 1 ), colour="Black")
  dMat <- rbind(dMat,c(Tempvar,names(GIlist)[i],seq(1.2,1.2+(length(GIlist)*0.1),0.1)[i], Datafinal$ST[bleh],Datafinal$Seqname[bleh],Datafinal$Dataset[bleh]))
}
colnames(dMat) <- c("x1","x2","Seq","y1","ST","Names","Dataset")

dMat$x1 <- as.numeric(dMat$x1)
dMat$x2 <- as.numeric(dMat$x2)
dMat$y1 <- as.numeric(dMat$y1)


#Change order of dMat so they're ordered by sequence type
Temp1 <- dMat$y1
dMat <- dMat[order(dMat$ST),]
dMat$y1 <- Temp1

#Prettier names for plotting purposes
dMat$dNames <- unlist(strsplit(unlist(strsplit(dMat$Names,"\\.\\d_"))[c(FALSE,TRUE)],"_genomic.fna.txt"))
# dSubjname <- unlist(strsplit(unlist(strsplit(basename(dSubject),"\\.\\d_"))[c(FALSE,TRUE)],"_genomic.fna"))

#Define PA7 serotype cluster positions and PA7 genome size
GIleft <- 2002773
GIright <- 2028627
dSubjsize <- 6588339

dMat$ST[dMat$ST == "9999"] <- "n/a"
dMat$ST <- as.factor(dMat$ST)
dMat$y1 <- rev(dMat$y1)

# Paper figure

#################################





library(ggplot2)
plotobj <- ggplot() +
  geom_rect(aes(xmin=GIleft,xmax=GIright,ymin=0.98,ymax=1.15+(0.1*length(GIlist))), color = "black", fill="black", alpha = 0.4) +
  geom_segment(data = dMat, aes(x=x1,xend=x2,y=y1,yend=y1), size = 0.8) +
  scale_alpha( guide="none")+
  geom_rect(data = dMat, aes(xmin = x1, xmax = x2, ymin = y1-0.05, ymax = y1+0.1-0.05 , fill = ST ),alpha = 0.4) +
  scale_size(guide="none") +
  annotate("text", x=mean(c(GIleft,GIright)),y=0.88, hjust="middle",label="O-antigen cluster",size=4) +
  # annotate("rect",xmin =1991539,xmax=1994304 ,ymin =1, ymax=8.705,alpha=0.5) + #GyrA
  # annotate("text",x =1992922, y = 0.95, hjust = "middle", label="gyrA", size = 4) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # axis.ticks.x = element_blank(),
  ) + 
  # theme_void() +
  coord_cartesian(xlim=c(min(1.8*1e06),(2.4*1e06) )  ,ylim=c(1,as.numeric(1.2+(length(GIlist)*0.1)))) +
  # coord_cartesian(xlim=c(min(dMat$x1),(max(dMat$x2)) )  ,ylim=c(1,as.numeric(1.2+(length(GIlist)*0.1)))) 
  # scale_x_continuous(name = "Genomic position in PA7",breaks=c(seq(min(dMat$x1),max(dMat$x2),200000))) 
  scale_x_continuous(name = "Genomic position in PA7",breaks=c(seq(1800000,2800000,200000)))

plotobj






# 
# plotobj2 <- ggplot(dMat) +
#   geom_rect(aes(xmin=GIleft,xmax=GIright,ymin=0.98,ymax=1.15+(0.1*length(GIlist)) ),alpha=0.03,fill="black") +
#   geom_segment(aes(x=x1,xend=x2,y=y1,yend=y1,colour=ST),size=1) +
#   scale_color_manual(values = yo, name = "ST") +
#   # scale_color_brewer(name = "Sequence type", palette = "Paired") +
#   # scale_size(guide="none") +
#   geom_text(aes(x=x1-1000,y=y1,label=dNames,hjust="right"),size=2) +
#   geom_segment(aes(x = 0, y = 1 , xend= dSubjsize, yend= 1), colour = "Black") + #PA7
#   annotate("text", x=mean(dMat$x1)-10000,y=1.07, label=dSubjname,size = 2) +
#   # geom_text(aes(x=mean(dMat$x1)-10000,y=1.07),label=dSubjname) +
#   # geom_rect(aes(xmin=GIleft,xmax=GIright,ymin=0.98,ymax=1.2+(0.1*length(GIlist)) ),alpha=0.03,fill="black") +
#   annotate("text", x=GIright+10000,y=1.07, label="O-antigen cluster",size = 2,hjust="left") +
#   # geom_text(aes(x=mean(c(GIleft,GIright)),y=1.07),label="O-antigen cluster") +
#   theme_minimal() +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#   ) + 
#   coord_cartesian(xlim=c(min(dMat$x1)-200000,(max(dMat$x2)) )  ,ylim=c(1,as.numeric(1.2+(length(GIlist)*0.1)))) +
#   # coord_cartesian(xlim=c(0,(max(dMat$x2)) )  ,ylim=c(1,as.numeric(1.2+(length(GIlist)*0.1)))) +
#   scale_x_continuous(name = "Genomic position") 
# 
# plotobj2
