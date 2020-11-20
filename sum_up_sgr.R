#Rscript
tryCatch({
  library(GenomicRanges)
},error= function(e){
  message('Package "GenomicRanges" and its dependencies are not installed or not in .libPaths\n ','Found libPath:\n',.libPaths(),'\n')
  stop(e, call.=FALSE)}
)


##   Load all sgr in folder and sum them up, outputs a rda file with 
args <- commandArgs(TRUE)

cat('\n\n\nThis script sums up coverage traces usally resulting from "bedtools genomecov"\n\n ')

if (length(args)!=2) {
  stop("Usage: Rscript sum_up_sgr path/to/sgr_folder output_file.rda", call.=FALSE)
} 
path2sgr<-args[1] 
output.file <- args[2]

tryCatch({
  setwd(path2sgr)
}, error= function(e){
  message("Path or Folder does not exist!")
  stop(e, call.=FALSE)
})

###### Functions

###############################################################
make_list_from_sgr2 <- function(sgr){
  # sgr V1=chr, V2=pos, V3= value
  chr<-unique(sgr[,1])
  # check if first nucleotide is 1 in 2nd column
  sgr_list<- vector("list", length(chr)) 
  names(sgr_list) <- chr
  for(i in 1:length(chr)){
    v1 <- sgr[,3][sgr[,1] == chr[i]]
    v2 <- sgr[,2][sgr[,1] == chr[i]]
    #print(min(v2))
    sgr_list[[i]]<- v1
  }
  return(sgr_list)
}

########################################################
##  run script: 


files <- system('ls *.sgr',intern = T)
cat(length(files),' Files found:\n\n', files)

#################################
###################################
n <- 0
i <- 1
try({
  Tr1 <- read.table(files[i],sep='\t',stringsAsFactors = F)
  Tr1_L <- make_list_from_sgr2(Tr1)
  n <- n+1
})
for(i in 2:length(files)){
  try({
    Tr2 <- read.table(files[i],sep='\t',stringsAsFactors = F)
    Tr2_L <- make_list_from_sgr2(Tr2)
    Tr1_L <- mapply('+', Tr1_L, Tr2_L)
    n <- n+1
  })
  
}

cat(n,' sgr files were loaded!\n')

f11 <- c(rep(1,20),0,rep(-1,20))

Tr1_L_f <- lapply(Tr1_L, filter,filter=(f11*0.1),circular=T)
Tr2_L_f <- lapply(Tr1_L_f, filter,filter=(f11*0.1),circular=T)
Tr_trace <- RleList(Tr1_L)
Tr_trace1st <- RleList(Tr1_L_f)  #first deviation
Tr_trace2nd <- RleList(Tr2_L_f)  # 2nd dev


save(Tr_trace, Tr_trace1st, Tr_trace2nd, file = output.file )

cat('3 Trace files were saved: raw, 1st deviation and 2nd deviation in: \n', output.file, '\n')


# cd /home/harald/Dropbox/shell_scripts
#  Rscript sum_up_sgr.R /media/harald/Disk_12/Monika/1260a7497b4b1b285e2ed83ab40de8bc/sgr/test test1.rda
 










