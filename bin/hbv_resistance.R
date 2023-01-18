#!/usr/bin/env Rscript

# TODO:
# Sekvensene må være riktig orientert først. Bruke MULTIFASTA fra ABACAS.
# Skal jeg aligne alle skaffolds hvis det er flere? Eller ta den lengste?


# Tanker:
# I referanse-sekvensen (alignet), finne posisjonen til første bokstav "A, T, G, C" (enten i "seq.reference" eller i alignment-objektet).
# Finne også siste bokstaven i refereanse-sekvensen (alignet).
# Bruke disse som koordinater til å subsette alignmentet. Dvs. dra ut bare den delen som aligner til referansen.
# Jeg vet at aminosyrene starter på første posisjon til referansen. Jeg kan deretter translatere derfra.

# En annen approach:
# Translatere referansen fra begynnelsen. Da vil jeg få aminosyreendringer, men ikke se dem på nucleotidnivå...

library(tidyverse)
library(jsonlite)
library(rjson)
#library(xml2)
library(seqinr)
#library(GenomicAlignments)
#library(msa)
#library(rphast)

rtDomain <- args[1]
blast    <- args[2]
              
# First create a vector of the reference sequence with each amino acid as an element
rtDomain_seq <- unlist(
  # split the sequence between every character and convert to a matrix
  str_split(
    read.fasta(rtDomain, seqtype = "AA", as.string = TRUE)$RT_domain[1],
    pattern = "",
    simplify = TRUE
    )
  )

# Read pairwise alignment from blastx
blast_json <- jsonlite::fromJSON("/home/jonr/Prosjekter/learning_nextflow/blast.json")
# Get the hsps and convert to tibble
query_aligned_seq <- unlist(
  str_split(
    # Extract the query hsp
    blast_json[["BlastOutput2"]][["report"]][["results"]][["bl2seq"]][[1]][["hits"]][[1]][["hsps"]][[1]][["qseq"]],
    pattern = "",
    simplify = TRUE
  )
)

# Find the differences between the aligned sequence and the reference
tmp <- as_tibble(rtDomain_seq) %>% 
  pivot_longer(everything()) %>% 
  mutate("position" = str_remove(name, "V")) %>% 
  rename("ref_aa" = "value") %>% 
  select(position, ref_aa)

tmp2 <- as_tibble(query_aligned_seq) %>% 
  pivot_longer(everything()) %>% 
  mutate("position" = str_remove(name, "V")) %>% 
  rename("sample_aa" = "value") %>% 
  select(position, sample_aa) %>% 
  # Convert everything to upper characters
  mutate(sample_aa = toupper(sample_aa))

# join the two and create mutations
df <- left_join(tmp, tmp2) %>% 
  # Create a column indicating when the reference and the sample is not identical
  mutate(mut = case_when(
    ref_aa != sample_aa ~ "mutated"
  )) %>% 
  filter(mut == "mutated") %>% 
  # Create the mutation names
  unite("Mutation", c(ref_aa, position, sample_aa), remove = FALSE, sep = "") %>% 
  select(-mut)


# Where are the differences
rtDomain_seq == query_aligned_seq
# Extract the differences
rtDomain_seq[rtDomain_seq != query_aligned_seq]
query_aligned_seq[query_aligned_seq == rtDomain_seq]
setdiff(rtDomain_seq, query_aligned_seq)
identical(rtDomain_seq, query_aligned_seq)

algorithm <- "Muscle"

fasta <- "/home/jonr/Prosjekter/learning_nextflow/HBV_results_rnaviral/6_abacas/HBV112022A1_C.abacas.MULTIFASTA.fa"
Ref_RT <- "/home/jonr/Prosjekter/learning_nextflow/Data/HBV_references/HBV_RT_domain.fa"

seq.list <- readDNAStringSet(c(fasta, Ref_RT))  

# Aligne sekvensen mot referansen 
seq.aln <- msa::msa(seq.list[c(1,2)], algorithm)
x <- DNAMultipleAlignment(seq.aln)
DNAStr <- as(x, "DNAStringSet")
#seq<-toupper(as.character(DNAStr[grep(samples.to.analyze[k],names(DNAStr))]))
# Get the sequence of the sequence
seq <- toupper(as.character(DNAStr[-grep(names(seq.list)[length(seq.list)],names(DNAStr))]))    

results <- as.data.frame(matrix(nrow = 1, ncol = 4))
colnames(results) <- c("Sample", "Deletions", "Frameshift", "Insertions")
results$Frameshift <- "NO"
results$Insertions <- "NO"


# Dra ut sekvens nr. 2 i lista fra alignmentet.
# as.character gjør denne om til en streng, ink. "dash"
# strsplit splitter mellom hver karakter
# unlist lager en vektor med hver karakter som element og navn på elementene blir sekvensnavn og posisjon i vektoren
seq.reference <- unlist(base::strsplit(as.character(DNAStr[grep(names(seq.list)[length(seq.list)],names(DNAStr))]),""))


library(DECIPHER)

# perform the alignment via the translations
# change NA to 1, 2 or 3 if the readingFrame is known
aligned <- AlignSeqs(seq.list)

# view the alignment in a browser (optional)
BrowseSeqs(aligned, highlight=0)

AAaligned <- AlignTranslation(seq.list,
                            readingFrame=1,
                            type="AAStringSet") # return AA or DNA sequences?

# view the alignment in a browser (optional)
BrowseSeqs(AAaligned, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned,
                file="<<REPLACE WITH PATH TO OUTPUT FASTA FILE>>")


### NACHO CODE

seq.reference <- unlist(base::strsplit(as.character(DNAStr[grep(names(seq.list)[length(seq.list)],names(DNAStr))]),""))
if(length(seq.reference[seq.reference=="-"])!=0){
  results$Insertions<-paste(as.numeric(which(seq.reference=="-")), collapse = " / ")
  
  ins.n<-length(seq.reference[seq.reference=="-"])
  ins.fs<-"YES"
  if(which(seq.reference=="-")==29904) ins.fs<-"NO" 
  if(length(which(as.numeric(which(seq.reference=="-")) %in% non.codding )) == length(which(seq.reference=="-"))) ins.fs<-"NO"
}else{
  ins.fs<-"NO"
  ins.n<-0
}


if(ins.n>0){
  fram.s.insertio<-ins.n%%3
  if(fram.s.insertio!=0){
    results$Frameshift<-"YES"
  }else{
    if(max(as.numeric(which(seq.reference=="-")))-min(as.numeric(which(seq.reference=="-")))> length(as.numeric(which(seq.reference=="-")))){
      results$Frameshift<-"YES"
    }
  }
}
if(ins.fs=="NO")results$Frameshift<-"NO"

out.df<-as.data.frame(matrix(data = NA, nrow = 36, ncol = 4))
colnames(out.df)<-c("Length", "Elements", "Positions","FS")
out.df$Length<-c(1:36)

for (i in 1:nrow(out.df)) {
  dummy<-unlist(base::strsplit(seq, paste("[A-Z]\\-{",i,"}[A-Z]",sep = "")))
  
  out.df$Elements[i]<-length(dummy)-1  
  if(length(dummy)>1){
    
    dummy<-dummy[-length(dummy)]
    characters<-as.numeric(nchar(dummy))
    characters[1]<-characters[1]+2
    if(length(dummy)>1){
      for(j in 2:(length(characters))){
        characters[j] <- characters[j-1]+characters[j] + 2 +i 
      }
    }
    
    if(ins.n>0){
      for (c in 1:length(characters)) {
        characters[c] <- characters[c]- length(which(as.numeric(which(seq.reference=="-"))<characters[c]))
      }
    }
    
    out.df$Positions[i]<-paste(characters,collapse = ";")
    if(length(which(characters %in% non.codding))==length(characters)){ 
      out.df$FS[i]<-"NO"
    }else{
      out.df$FS[i]<-"YES"  
    }
  }
  
  
}

out.df$To.out<-paste(out.df$Length, "[", out.df$Positions,"]", sep = "")
out.df$To.out[out.df$Elements==0]<-NA

results$Deletions<-paste(na.omit(out.df$To.out), collapse = " / ")
results$Sample<-samples.to.analyze
deletion.lengh<-out.df$Length[out.df$Elements!=0]%%3
if(length(which(deletion.lengh>0))>0 & length(which(out.df$FS=="YES"))>0){
  results$Frameshift<-"YES"
}

date <- gsub("-","",Sys.Date())

# Skip this - instead use the final.results object
#write.csv(final.results, paste(results.folder, date, "DeletionFinderResults.csv",sep = ""), row.names = FALSE)

####################################
### Deletion Finder end ###
####################################



genes <- read.csv(genelist)
non.codding <- c(1:29903)
for (i in 1:nrow(genes)) {
  non.codding<-non.codding[-which(non.codding %in% c(genes$start[i]:genes$end[i]))]
}


#Cleaning
deletion_results <- results

deletion_results$Frameshift[which(deletion_results$Deletions=="1[28271] / 3[21991] / 6[21765] / 9[11288]" & deletion_results$Insertions=="NO")]<-"NO"

#Check insertions and deletion region of genes
positions.to.test<-list()
for (i in 1:nrow(deletion_results)) {
  if(deletion_results$Insertions[i]=="NO" & deletion_results$Frameshift[i]=="YES"){
    dummy<-unlist(base::strsplit(deletion_results$Deletions[i],"/"))
    del.check<-FALSE
    for (j in 1:length(dummy)) {
      size<-as.numeric(gsub("\\[.*", "",dummy[j]))
      
      if(size%%3 !=0){
        positions.del<-gsub("\\].*","",gsub(".*\\[","",dummy[j]))
        positions.del<-as.numeric(unlist(base::strsplit(positions.del,";")))
        if(length(positions.del)==length(positions.del[which(positions.del %in% non.codding)]) & deletion_results$Insertions[i]=="NO"){
          if(!del.check) deletion_results$Frameshift[i]<-"NO"
          rm(positions.del)
        }else{
          if(length(positions.del[which(positions.del %in% non.codding)])>0) positions.del<-positions.del[-which(positions.del %in% non.codding)]
          if(length(positions.to.test)==0) positions.to.test[[1]]<-positions.del
          if(length(positions.to.test)>0)positions.to.test[[length(positions.to.test)]]<-c(positions.to.test[[length(positions.to.test)]],positions.del)
          deletion_results$Frameshift[i]<-"YES"
          del.check<-TRUE
        }
      }
      
      
    }
    if(deletion_results$Insertions[i]=="NO" & deletion_results$Frameshift[i]=="YES"){
      names(positions.to.test)[length(positions.to.test)]<-deletion_results$Sample[i]
    }  
  }
  
}

deletion_results<-deletion_results[order(deletion_results$Frameshift, decreasing = TRUE),]
# Skip this as well
#write_xlsx(deletion_results[,c(1:4)],paste(results.folder,"FrameShift_", gsub("\\.fa.*","",gsub(".*/","", total.fasta)),".xlsx",sep = ""))


# FrameshiftDB ------------------------------------------------------------

indels<-read.csv(database)
df <- deletion_results

# If the sample has NO frameshifts then it is ready
if(length(which(df$Frameshift=="NO"))>0) {
  df.ready <- df[which(df$Frameshift=="NO"),]
  df.ready$Ready<-"YES"
  df.ready$Comments<-"No frameshifts detected"
}

if(length(which(df$Frameshift=="YES"))>0){
  
  df<-df[which(df$Frameshift=="YES"),]
  df$Ready<-"NO"
  df$Comments<-NA
  
  indels<-indels[which(indels$Status=="Confirmed Fastq"),]
  
  insertion.list<-gsub(".*: ","",indels$ID[grep("Insertion",indels$ID)])
  deletion.list<-gsub(".*: ","",indels$ID[grep("Deletion",indels$ID)])
  
  for (i in 1:nrow(df)) {
    flag<-NA
    dummy.ins<-gsub(" ","",unlist(base::strsplit(df$Insertions[i],"/")))
    dummy.dels<-gsub(" ","",unlist(base::strsplit(df$Deletions[i],"/")))
    
    if(length(dummy.ins[which(dummy.ins %in% insertion.list)])!=length(dummy.ins) & dummy.ins[1]!="NO"){
      
      if(length(dummy.ins)>1){ 
        df$Comments[i]<-paste("Unknown insertion/s detected at", paste(dummy.ins[-which(dummy.ins %in% insertion.list)], collapse = ";"))
        if(length(which(dummy.ins %in% insertion.list))==0){
          df$Comments[i]<-paste("Unknown insertion/s detected at", paste(dummy.ins,collapse = ","))
        }
        
      }else{
        df$Comments[i]<-paste("Unknown insertion/s detected at", dummy.ins)
      }
      flag<-"InsKO"
    }else{
      flag<-"InsOK"
    }
    
    to.clean<-dummy.dels[grep(";", dummy.dels)]
    
    if(length(to.clean)>0){
      dummy.dels<-dummy.dels[-grep(";", dummy.dels)]
      for (j in 1:length(to.clean)) {
        dummy.dels2<-unlist(base::strsplit(to.clean[j],";"))
        dummy.dels2[-1]<- paste(gsub("\\[.*","[",dummy.dels2[1]), dummy.dels2[-1],sep = "")
        dummy.dels2<-paste(dummy.dels2,"]",sep = "")
        dummy.dels2<-gsub("]]","]",dummy.dels2)
        dummy.dels<-c(dummy.dels2, dummy.dels)
      }
    }
    dummy.dels<-dummy.dels[which(as.numeric(gsub("\\[.*","",dummy.dels))%%3 !=0 )]
    
    if(length(dummy.dels)==0 & is.na(df$Comments[i]) ){
      #No deletions FS and all Insertions are OK
      df$Ready[i]<-"YES"
      df$Comments[i]<-"All frameshifts are OK"
    }
    
    #INS OK / DELS Ok
    if(length(dummy.dels)>0){
      if(length(dummy.dels[which(dummy.dels %in% deletion.list)])==length(dummy.dels) & flag=="InsOK" ){
        df$Ready[i]<-"YES"
        df$Comments[i]<-"All frameshifts are OK"
      }
      
      #INS OK /DELS KO
      if(length(dummy.dels[which(dummy.dels %in% deletion.list)])!=length(dummy.dels) & flag=="InsOK" ){
        df$Ready[i]<-"NO"
        if(length(dummy.dels)>1){ 
          if(length(which(dummy.dels %in% deletion.list))>0){
            df$Comments[i]<-paste("Unknown deletions/s detected at", paste(dummy.dels[-which(dummy.dels %in% deletion.list)], collapse = ";"))
          }else{
            df$Comments[i]<-paste("Unknown deletions/s detected at", paste(dummy.dels,collapse = ","))
          }
        }else{
          df$Comments[i]<-paste("Unknown deletions/s detected at", dummy.dels)
        }
        flag<-"InsOK_DelKO"
      }
      
      #INS KO /DELS KO
      if(length(dummy.dels[which(dummy.dels %in% deletion.list)])!=length(dummy.dels) & flag=="InsKO" ){
        df$Ready[i]<-"NO"
        if(length(dummy.dels)>1){ 
          df$Comments[i]<-paste(df$Comments[i], "&", paste("Unknown deletions/s detected at", paste(dummy.dels[-which(dummy.dels %in% deletion.list)], collapse = ";")))
        }else{
          df$Comments[i]<-paste(df$Comments[i], "&", paste("Unknown deletions/s detected at", dummy.dels))
        }
      }
      
      
      
      
      
