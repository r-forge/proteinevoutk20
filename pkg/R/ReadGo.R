##the order of GO names in Rokas's data
gonames <- read.csv("~/proteinevoutk20/pkg/Data/GoNames.csv",header=FALSE)
##the groups according to GO, functions of genes
filename <- "~/proteinevoutk20/pkg/Result/slimMapperResult-Function.xlsx"
#filename <- "~/proteinevoutk20/pkg/Result/slimMapperResult-YeastComponent.xlsx"
GoGroup <- read.xlsx(filename,sheetIndex=1,header=TRUE,colIndex=c(2,6))
##GoGroupList
GGList <- vector(mode="list",length=dim(GoGroup)[1])
## the following convert GO names in each group to the order of genes in Rokas's data
for(i in 1:dim(GoGroup)[1]){
  go <- unlist(strsplit(toString(GoGroup[i,2]),split=","))
  GGList[[i]]$GOterm <- toString(GoGroup[i,1]) #GO term(depends on function, component or others)
  GGList[[i]]$genes <- go # GO names in this group
  GGList[[i]]$order <- which(gonames[,2] %in% go) # order of genes in this group
}