
########Costumize the whole script by uncommenting according to the instructions######## 


########Install microhaplot (uncomment for installing the package and moving the files the first time you run the script) ######## 

#install.packages("devtools")
#devtools::install_github("ngthomas/microhaplot", build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
#microhaplot::mvShinyHaplot("~/Shiny") 




######## Run the perl script to get raw haplotype calls ######## 


#Costumize RUN_NAME and directory path in case the folder is not in Desktop

run.label <- "S5_VISAGE_ET_MICROHAPLOTYPER"
sam.path <- "~/Desktop/S5_VISAGE_ET_MICROHAPLOTYPER"
label.path <- "~/Desktop/S5_VISAGE_ET_MICROHAPLOTYPER/label.txt"
vcf.path <- "~/Desktop/S5_VISAGE_ET_MICROHAPLOTYPER/CALLS.vcf"
app.path <- "~/Shiny/microhaplot"


#######Customize thresholds to apply for filtering######
### Minimum allele frequency threshold (0 to 1)
min.al.freq=0.02    
### Minimum coverage per alelle
min.cov=15

#######Uncomment to check your microhaplot version######

packageVersion("microhaplot")

library (microhaplot)

#######Uncomment to run microhaplot version 0.1.0 ######

haplo.read.tbl <- runHaplot(run.label = run.label,
                            sam.path=sam.path,
                            label.path=label.path,
                            vcf.path=vcf.path,
                            app.path=app.path)


#######Uncomment to run microhaplot version 1.0.0.9000 ######

haplo.read.tbl <- prepHaplotFiles(run.label = run.label,
                            sam.path=sam.path,
                            label.path=label.path,
                            vcf.path=vcf.path,
                            app.path=app.path)




######## Save raw data #######

library(tidyverse)
locinames <- c("1pA", "MH01", "1pD", "MH03", "MH04", "3pB", "3qC", "4qD", "7pB", "8pA", "8pB", "9pA", "10pB", "MH11", "12qB", "15qD", "MH18", "MH20", "MH21", "MH22", "22qB")
haplo.read.tbl$locus <- factor(haplo.read.tbl$locus, levels = locinames)
haplo.read.tbl <- haplo.read.tbl %>% arrange(id,locus,desc(depth))
write.csv(haplo.read.tbl, file="haplo.read.tbl.txt" )


######## Start from raw data ######## 
library(dplyr)
haplo.read.tbl<-read.csv("haplo.read.tbl.txt")
locinames <- c("1pA", "MH01", "1pD", "MH03", "MH04", "3pB", "3qC", "4qD", "7pB", "8pA", "8pB", "9pA", "10pB", "MH11", "12qB", "15qD", "MH18", "MH20", "MH21", "MH22", "22qB")
haplo.read.tbl$locus <- factor(haplo.read.tbl$locus, levels = locinames)
haplo.read.tbl <- haplo.read.tbl %>% arrange(id,locus,desc(depth))


######## Filter calls ######## 

library(dplyr)

FilteringTable <- haplo.read.tbl %>% filter(!grepl("N",haplo)) %>% filter(!grepl("X",haplo))

preTotaldepth<-FilteringTable %>% 
  group_by(locus,id) %>% 
  summarise(preTotaldepth = sum(depth))

FilteringTable <- merge(FilteringTable, preTotaldepth, all = TRUE)
FilteringTable<- FilteringTable %>% mutate(preAl.freq=depth/preTotaldepth)
FilteredTable= FilteringTable %>% filter(depth>=min.cov & preAl.freq>=min.al.freq)
FilteredTable <- FilteredTable %>% select(-c(preAl.freq,preTotaldepth))

Totaldepth<-FilteredTable %>% 
  group_by(locus,id) %>% 
  summarise(Totaldepth = sum(depth))

FilteredTable <- merge(FilteredTable, Totaldepth, all = TRUE)

FilteredTable<- FilteredTable %>% mutate(Al.freq=depth/Totaldepth)

MisperformingMarkers <- c()
`%not_in%` <- purrr::negate(`%in%`)
FilteredTable <- droplevels(FilteredTable %>% filter (locus %not_in% MisperformingMarkers))




######## Split data and get profile plots######## 
data_split <- split(FilteredTable, FilteredTable$id)
Samples <- as.character(paste("Sample",unique(FilteredTable$id),sep = "_"))

for (i in 1:length(data_split)) {
  assign(Samples[i], data_split[[i]])}

profile.plot <- function(data, title)  {
   ggplot(data,aes(x=haplo,y=depth))+
    geom_bar(stat="identity")+
    facet_wrap(~locus,scales="free")+
    theme(axis.text.x = element_text(hjust = 0.5, vjust= 0, size= 6), axis.text.y=element_text(size=4))+
    ggtitle(title)+
    theme(plot.title = element_text(size=12, face="bold",  hjust=0))+
    theme(axis.title = element_text(face="bold", size= 12))+
    ylab("Coverage") + xlab("MH-allele")+
    theme(strip.text.x = element_text(size = 7, margin = margin(.05, 0, .05, 0, "cm"),face="bold"))}

for(i in Samples){
  myplot<-(profile.plot(get(i), i))
  ggsave(myplot,filename=paste(i,".pdf",sep=""), width = 14, height = 8.5, units = "in")} 


#####Get genotype tables#####

library(tidyr)
library(data.table)

locinames <- c("1pA", "MH01", "1pD", "MH03", "MH04", "3pB", "3qC", "4qD", "7pB", "8pA", "8pB", "9pA", "10pB", "MH11", "12qB", "15qD", "MH18", "MH20", "MH21", "MH22", "22qB")
for(i in Samples){
  GenotypingTable<-select(get(i),locus,haplo,depth)
  GenotypingTable$locus <- factor(GenotypingTable$locus, levels = locinames)
  OrderGenotypeTable<-GenotypingTable %>% group_by(locus) %>% mutate(haplo.number = order(order(locus, depth, decreasing=TRUE)))
  GenotypeTable<-dcast(setDT(OrderGenotypeTable), locus~haplo.number, value.var=c('haplo', 'depth'))
  GenotypeTable2<-if("haplo_3" %in% colnames(GenotypeTable)){
    GenotypeTable%>%mutate(genotype=case_when(
        !is.na(haplo_3)  ~ paste0(haplo_1,"/",haplo_2,"*"),
        is.na(haplo_2)  ~ paste0(haplo_1,"/",haplo_1),
        !is.na(haplo_2) ~ paste0(haplo_1,"/",haplo_2)))
  }else{
    GenotypeTable%>%mutate(genotype=case_when(
      is.na(haplo_2)  ~ paste0(haplo_1,"/",haplo_1),
      !is.na(haplo_2) ~ paste0(haplo_1,"/",haplo_2)))}
  write.csv(GenotypeTable2, file=paste(i,"_genotypes.txt",sep=""))} 
