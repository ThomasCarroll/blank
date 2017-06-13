library(GenomicFeatures)
library(rtracklayer)

ChIPQCAnnotationFromGFF3 <- function(GFF3,GeneAnnotation="Custom"){
  txdbFromGFF3 <- makeTxDbFromGFF(GFF3)
  return(ChIPQCAnnotationFromTXDB(txdbFromGFF3,GeneAnnotation))
}
ChIPQCAnnotationFromTXDB <- function(txdb,GeneAnnotation="Custom",AllChr=NULL){
  All5utrs <- reduce(unique(unlist(fiveUTRsByTranscript(txdb))))
  All3utrs <- reduce(unique(unlist(threeUTRsByTranscript(txdb))))
  Allcds <- reduce(unique(unlist(cdsBy(txdb,"tx"))))
  Allintrons <- reduce(unique(unlist(intronsByTranscript(txdb))))
  Alltranscripts <- reduce(unique(unlist(transcripts(txdb))))
  
  posAllTranscripts <- Alltranscripts[strand(Alltranscripts) == "+"]
  posAllTranscripts <- posAllTranscripts[!(start(posAllTranscripts)-20000 < 0)]
  negAllTranscripts <- Alltranscripts[strand(Alltranscripts) == "-"]
  chrLimits <- seqlengths(negAllTranscripts)[as.character(seqnames(negAllTranscripts))]
  if(!any(is.na(chrLimits))){
    negAllTranscripts <- negAllTranscripts[!(end(negAllTranscripts)+20000 > chrLimits)] 
  }
  Alltranscripts <- c(posAllTranscripts,negAllTranscripts)
  Promoters500 <-  reduce(flank(Alltranscripts,500))    
  Promoters2000to500 <-  reduce(flank(Promoters500,1500))
  LongPromoter20000to2000  <- reduce(flank(Promoters2000to500,18000))
  if(!is.null(AllChr)){
    All5utrs <- GetGRanges(All5utrs,AllChr=AllChr)
    All3utrs <- GetGRanges(All3utrs,AllChr=AllChr)
    Allcds <- GetGRanges(Allcds,AllChr=AllChr)
    Allintrons <- GetGRanges(Allintrons,AllChr=AllChr)
    Alltranscripts <- GetGRanges(Alltranscripts,AllChr=AllChr)
    Promoters500 <- GetGRanges(Promoters500,AllChr=AllChr)
    Promoters2000to500 <-  GetGRanges(Promoters2000to500,AllChr=AllChr)
    LongPromoter20000to2000  <- GetGRanges(LongPromoter20000to2000,AllChr=AllChr)        
  }
  return(list(version=GeneAnnotation,LongPromoter20000to2000=LongPromoter20000to2000,
              Promoters2000to500=Promoters2000to500,Promoters500=Promoters500,
              All5utrs=All5utrs,Alltranscripts=Alltranscripts,Allcds=Allcds,
              Allintrons=Allintrons,All3utrs=All3utrs))
}

