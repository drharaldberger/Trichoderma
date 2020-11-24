#Rscript
# Load the rda file and calculate transcription start point (TSP) 

tryCatch({
  library(GenomicRanges)
  library(stringr)
},error= function(e){
  message('Package "GenomicRanges" and "stringr" or their dependencies are not installed or not in .libPaths\n ','Found libPath:\n',paste(.libPaths(),collapse = '\n'),'\n')
  stop(e, call.=FALSE)}
)

args <- commandArgs(TRUE)

cat('\n\n\nThis script estimates transcription start point (TSP) \nbased on transcript traces;\n it requires a gff file with exon features\n\n ')

if (length(args)!=5) {
  stop("Usage: Rscript TSP_from_transcript input.rda input.gff e.g.:Parent=mRNA_ genome.txt output_file.txt transcript_threshold")
} 
input_file <-args[1]  #rda file
gff_file <- args[2]
gene_ID <- args[3]  ## what needs to be removed
#genome_sizes <- args[4]
output.file <- args[5]
tr_threshold <- as.numeric(args[4])
## Rscript /home/harald/Dropbox/shell_scripts/TSP_from_transcript.R test1.rda /home/harald/Dropbox/Monika/From_Wolfgang_new_Trichoderma_anno/QM6a_gene_annotation/T_reesei_nG_2.gff3 Parent=mRNA_ /home/harald/Dropbox/Monika/From_Wolfgang_new_Trichoderma_anno/QM6a_gene_annotation/Chromsizes.txt test2.bed 0
#input_file <- '/media/harald/Disk_12/Monika/1260a7497b4b1b285e2ed83ab40de8bc/sgr/rev/rev_all.rda'  #rda file
#input_file <- '/media/harald/Disk_12/Monika/1260a7497b4b1b285e2ed83ab40de8bc/sgr/forw/forw_all.rda'  #rda file

#gff_file <- '/home/harald/Dropbox/Monika/From_Wolfgang_new_Trichoderma_anno/QM6a_gene_annotation/T_reesei_nG_2.gff3'
#gene_ID <- 'Parent=mRNA_'  ## what needs to be removed
#genome_sizes <- '/home/harald/Dropbox/Monika/From_Wolfgang_new_Trichoderma_anno/QM6a_gene_annotation/Chromsizes.txt'
#output.file <- 'test2.bed'
#tr_threshold <- as.numeric('0')
#setwd('/media/harald/Disk_12/Monika/1260a7497b4b1b285e2ed83ab40de8bc/sgr/DD_f')

tryCatch({
  load(input_file, verbose = T)  # Tr_trace, Tr_trace1st, Tr_trace2nd
}, error= function(e){
  message(".rda file does not exist!")
  stop(e)
})
cat('.rda file loaded\n\n')
cat('\n\n')
tryCatch({
  gff <- read.table(gff_file,sep='\t',stringsAsFactors = F)
  
}, error= function(e){
  message(".gff file does not exist!")
  stop(e)
})
cat('.gff file loaded\n\n')


#gff:
# unitig_0_consensus	maker	exon	7295	9560	.	+	.	Parent=mRNA_TrA0001W

gff2 <- gff[gff$V3=='exon',]
exons <- sub(gene_ID,'',gff2$V9)

gff_ir <- IRanges(start = gff2$V4, end = gff2$V5, names = exons)
gff_gr <- GRanges(seqnames = gff2$V1, ranges= gff_ir, strand = gff2$V7)

##  filter gff for sequences not in the traces
#sel_gff <- as.vector(seqnames(gff_gr)) %in% names(Tr_trace) 
#gff_gr <- gff_gr[sel_gff]
## get the first and the last exon per gene
## make genome sizes from Tr_trace

chr_len <- sapply(Tr_trace,length)
genome.sizesGR <- GRanges(seqnames = names(Tr_trace), ranges = IRanges(start = 1, end = chr_len))

### allow only genes wich are covered
sel_1 <- gff_gr %within% genome.sizesGR
gff_gr <- gff_gr[sel_1]

genes <- unique(names(gff_gr))
cat('Calculate exon transcript levels...\n\n')
min_Gr <- GRanges()
max_Gr <- GRanges()

for(i in 1:length(genes)){
  test1 <- gff_gr[names(gff_gr)==genes[i]]
  min1 <- test1[min(start(test1)) == start(test1)]
  min_Gr <- c(min_Gr,min1)
  max1 <- test1[max(start(test1)) == start(test1)]
  max_Gr <- c(max_Gr,max1)
  
}
cat('5prime and 3prime exons identified.\n\n')

##count zeros in exons

count_zeros <- function(Rle_in){
  sum(Rle_in == 0)
}
###
min_Rle <- Tr_trace[min_Gr]
mean_min_Rle <- sapply(min_Rle,median)
count_Z_min_Rle <- sapply(min_Rle,count_zeros)
length_min_Rle <- sapply(min_Rle,length)
max_Rle <- Tr_trace[max_Gr]
mean_max_Rle <- sapply(max_Rle,median)
count_Z_max_Rle <- sapply(max_Rle,count_zeros)
length_max_Rle <- sapply(max_Rle,length)
min_max_df <- data.frame(min= mean_min_Rle, minZ=count_Z_min_Rle, min_len=length_min_Rle, max= mean_max_Rle,maxZ=count_Z_max_Rle, max_len=length_max_Rle)

cat('\nTranscription threshold provided: ',tr_threshold,'\n')
## find threshold
if(tr_threshold==0){
  exZmin <- min_max_df$minZ == 0 & min_max_df$min_len > 500
  exZmax <- min_max_df$maxZ == 0 & min_max_df$max_len > 500
   
  q1 <- quantile(c(min_max_df$max[exZmax],min_max_df$min[exZmin]),probs = seq(0, 1, 0.1))
  #barplot(q1,ylim = c(0,50))
  q11 <- q1[2]
  
  if(q11==0){
    for(i2 in 3:length(q1)){
      q11 <- q1[i2]
      if(q11!=0){
        cat('Transcript threshold value set to:',q11,'\n according to',names(q11),'quantile\n\n')
      }
    }
  }else{
    cat('Transcript threshold value set to: ',q11,'\n according to 10% quantile\n\n')
  }
  tr_threshold <- q11
}



############

detect_TSP <- function(Tr_trace, Tr_trace1st ,Tr_trace2nd , exon5prime, exon3prime, exon_tr, thresh=tr_threshold, genome.sizesGR){
  # length(exon5prime) == length(exon3prime) == nrow(exon_tr)
  # outputs
  #TSP<- c(); TSP_p <- c(); TSP_desc <- c(); genes1 <- c();
  strand2 <- genes1 <- ATG <- TSP_desc <- TSP_p <- TSP <- rep(NA,length(exon5prime))
  for(i in 1:length(exon5prime)){
    strand1 <- as.character(strand(exon5prime[i]))
    strand2[i] <- strand1
    gene1 <- names(exon5prime[i])
    cat(' ',gene1)
    if(strand1 =='+'){
      ## Forward
      exon_tr1 <- exon_tr[i,1]
      if(exon_tr1 < thresh){
        # skip if not enough transcript available
        TSP[i] <- start(exon5prime[i])      ###################output if no transcript is found
        ATG[i] <- start(exon5prime[i])
        TSP_p[i] <- 0
        TSP_desc[i] <- 'no Transcript'
        genes1[i] <- gene1
        next
      }
      start1 <- start(exon5prime[i])
      out_slopes <- c()
      out_pos_slopes <- c()
      for(i2 in 1:50){
        #go down in 40 bp steps, max 50 steps
        steps <- 40
        
        gr_sel <- GRanges(seqnames = seqnames(exon5prime[i]), ranges = IRanges(start= start1-steps,width=steps+1))
        
        ## check if still in genome
        if(gr_sel %within% genome.sizesGR){
          if(gr_sel %over% exon5prime[-i] | gr_sel %over% exon3prime[-i]){
            break
          }
          v1 <- suppressWarnings(as.numeric(unlist(Tr_trace2nd[gr_sel])))
          dec1 <-descent_through(v1)
          out_pos_slopes <- c(which(dec1)+(start1-steps), out_pos_slopes)
          ## the slope (1st dev give the probability)
          v2 <- suppressWarnings(as.numeric(unlist(Tr_trace1st[gr_sel])))
          out_slopes <- c(v2[dec1], out_slopes)   #append outslopes
          v3 <- suppressWarnings(as.numeric(unlist(Tr_trace[gr_sel])))
          
          if((mean(v3)< exon_tr1 / 5) | (mean(v3)> exon_tr1 * 2)  ){
            break
          }
          
          # which(genes=='TrA1039W')
        }else{
          break
        }
        start1 <- start1 - steps ## going back 5'
      } 
      if(length(out_slopes)==0){
        TSP[i] <- start(exon5prime[i])     #########  outout if no slopes are detected
        ATG[i] <- start(exon5prime[i])
        TSP_p[i] <- 0
        TSP_desc[i] <- 'no TSP detected'
        genes1[i] <- gene1
      }else{
        # take the most 5' one & take only if the slope is at least a quater of the max slope in the 5' region
        sel1 <- out_slopes > (max(out_slopes)/4) 
        out_pos_slopes <- out_pos_slopes[sel1]
        out_slopes <- out_slopes[sel1]
        if(length(out_slopes)==0){
          TSP[i] <- start(exon5prime[i])     #########  outout if no slopes are detected
          ATG[i] <- start(exon5prime[i])
          TSP_p[i] <- 0
          TSP_desc[i] <- 'no TSP detected'
          genes1[i] <- gene1
        }else{
          TSP[i] <- out_pos_slopes[1]
          ATG[i] <- start(exon5prime[i])
          TSP_p[i] <- (out_slopes[1]/(exon_tr1*2 + thresh))*1000 #relation of slope vs exon
          TSP_desc[i] <- 'TSP detected'
          genes1[i] <- gene1
        }
        
      }
      
    }else{
      ## reverse
      
      exon_tr1 <- exon_tr[i,4]   ##changed
      if(exon_tr1 < thresh){
        # skip if not enough transcript available
        TSP[i] <-  end(exon3prime[i])      ##########  output if no transcript
        ATG[i] <-  end(exon3prime[i])
        TSP_p[i] <- 0
        TSP_desc[i] <- 'no Transcript'
        genes1[i] <- gene1
        next
      }
      start1 <- end(exon3prime[i])   ##changed
      out_slopes <- c()
      out_pos_slopes <- c()
      for(i2 in 1:50){
        #go down in 20 bp steps, max 50 steps
        steps <- 40
        
        gr_sel <- GRanges(seqnames = seqnames(exon3prime[i]), ranges = IRanges(start= start1, width=steps+1))
        
        ## check if still in genome
        if(gr_sel %within% genome.sizesGR){
          if(gr_sel %over% exon5prime[-i] | gr_sel %over% exon3prime[-i]){
            break
          }
          v1 <- suppressWarnings(as.numeric(unlist(Tr_trace2nd[gr_sel])))
          as1 <-ascent_through(v1)
          out_pos_slopes <- c( out_pos_slopes, which(as1)+(start1))
          ## the slope (1st dev give the probability)
          v2 <- suppressWarnings(as.numeric(unlist(Tr_trace1st[gr_sel])))
          out_slopes <- c( out_slopes, v2[as1])   #append outslopes
          v3 <- suppressWarnings(as.numeric(unlist(Tr_trace[gr_sel])))
          
          if((mean(v3) < exon_tr1 / 5)  | (mean(v3)> exon_tr1 * 2) ){
            break
          }
        }else{
          break
        }
        start1 <- start1 + steps ## going forward 3'
      } 
      if(length(out_slopes)==0){
        TSP[i] <-  end(exon3prime[i])    ####  output if no slppes are detected
        ATG[i] <-  end(exon3prime[i])
        TSP_p[i] <- 0
        TSP_desc[i] <- 'no TSP detected'
        genes1[i] <- gene1
      }else{
        # take the most 5' one & take only if the slope is at least a quater of the max slope in the 5' region
        sel1 <- out_slopes < (min(out_slopes)/4) 
        out_pos_slopes <- out_pos_slopes[sel1]
        out_slopes <- out_slopes[sel1]
        if(length(out_slopes)==0){
          TSP[i] <-  end(exon3prime[i])    ####  output if no slppes are detected
          ATG[i] <-  end(exon3prime[i])
          TSP_p[i] <- 0
          TSP_desc[i] <- 'no TSP detected'
          genes1[i] <- gene1
        }else{
          TSP[i] <- out_pos_slopes[length(out_pos_slopes)]
          ATG[i] <-  end(exon3prime[i])
          TSP_p[i] <- abs(out_slopes[length(out_slopes)]/(exon_tr1*2 + thresh))*1000 #relation of slope vs exon
          TSP_desc[i] <- 'TSP detected'
          genes1[i] <- gene1
        }
        
      }
      
    }
  }
  return(data.frame(genes=genes1, TSP, ATG, strand=strand2, TSP_score=TSP_p, TSP_desc, stringsAsFactors = F))
}



###############################################
descent_through<- function(x,thresh=0){
  above<- x > thresh
  intersect.points<-which(diff(above) < 0)  #gives -1 or 1 for points before the intersect
  x.slopes<-x[intersect.points+1]-x[intersect.points]
  
  x.points<-((thresh-x[intersect.points])/x.slopes)+intersect.points
  x.point<-round(x.points)
  out1<-rep(F,length(x))
  out1[x.point]<-T  ## make true false vector
  return(out1)
}
################################################
###############################################
ascent_through<- function(x,thresh=0){
  above<- x > thresh
  intersect.points<-which(diff(above) > 0)  #gives -1 or 1 for points before the intersect
  x.slopes<-x[intersect.points+1]-x[intersect.points]
  
  x.points<-((thresh-x[intersect.points])/x.slopes)+intersect.points
  x.point<-round(x.points)
  out1<-rep(F,length(x))
  out1[x.point]<-T  ## make true false vector
  return(out1)
}
################################################



cat('TSP detection started....\n\n')
cat('Genes analysed:\n')
TSP.df.test <- detect_TSP(Tr_trace=Tr_trace, Tr_trace1st=Tr_trace1st ,Tr_trace2nd=Tr_trace2nd, exon5prime=min_Gr, exon3prime=max_Gr, exon_tr=min_max_df, genome.sizesGR=genome.sizesGR)
cat('\n\nDone!!!\n')
m1<- apply(TSP.df.test[,2:3],1,sort)
df.TSP.test <- data.frame(genes=TSP.df.test$genes,t(m1),TSP.df.test$strand,round(TSP.df.test$TSP_score), stringsAsFactors = F)


df.genes <- data.frame(chr=gff2$V1,genes=exons,stringsAsFactors = F)
df.genesU <- unique(df.genes)

df.TSP.merged <- merge(df.genesU,df.TSP.test,by='genes')
df.TSP.merged <- format(df.TSP.merged, scientific=F, trim = T)
write.table(df.TSP.merged[,c(2:4,1,6,5)],output.file, sep = '\t', row.names = F,col.names = F,quote = F)
cat('\n\nWrite .bed file:\n',output.file,'\n\n')

##################################################################################################################
