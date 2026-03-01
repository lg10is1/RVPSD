
library("tidyverse")
virus_name="$name" #rename
wd=paste("$path_to_dir/",virus_name,"/",sep="")
wd1=paste("path_to_dir/result_best_pdb/",virus_name,"/",sep="")
wd2=paste("path_to_dir/fasta_file_dir/",virus_name,"/",sep="") 
setwd(wd)
dir.create(wd1)

d1=list.files(wd)

  
f=data.frame(c("Failed_proteins"),"length", "Fasta_Path")
plddt.file1=data.frame(c("Protein_accession"),"plddt")
for (n in 1:length(d1)) {

  d2=list.files(str_c(wd,d1[n]))

  if ("ranked_0.pdb" %in%d2){

   n12=d1[n]


    file.copy(str_c(wd,d1[n],"/ranked_0.pdb"),str_c(wd1,n12,".pdb"))
    plddt1=read_lines(str_c(wd,d1[n],"/ranking_debug.json"))[3:7]
    plddt1.1=gsub(",","",plddt1)
    plddt1.2=gsub("        \"model_1_pred_0\": ","",plddt1.1)
    plddt1.3=gsub("        \"model_2_pred_0\": ","",plddt1.2)
    plddt1.4=gsub("        \"model_3_pred_0\": ","",plddt1.3)
    plddt1.5=gsub("        \"model_4_pred_0\": ","",plddt1.4)
    plddt1.6=gsub("        \"model_5_pred_0\": ","",plddt1.5)
    plddt2=max(as.numeric(plddt1.6)) 

    plddt.file=data.frame(n12,plddt2)
    colnames(plddt.file)=c("Protein_accession","plddt")

    colnames(plddt.file1)=c("Protein_accession","plddt")

    plddt.file1=rbind(plddt.file1,plddt.file)
	
  }else{
    fa=read_lines(str_c(wd2,d1[n],".fasta"))
	fasta_path = str_c(wd2,d1[n],".fasta")
	f.1=data.frame(d1[n],sum(nchar(fa[2:length(fa)])),fasta_path)
	    colnames(f)=c("Failed_proteins","length","Fasta_Path")

    colnames(f.1)=c("Failed_proteins","length","Fasta_Path")
	f=rbind(f,f.1)
]
  }
}
f1=as.data.frame(f)
    colnames(f1)=c("Failed_proteins","length","Fasta_Path")
    f1=na.omit(f1)
	f1=f1[-1,]
  write.csv(f1,str_c(wd1,"failed_proteins.csv"),row.names = F)
  
  plddt.file1=na.omit(plddt.file1)
  plddt.file1=plddt.file1[-1,]
  write.csv(plddt.file1,str_c(wd1,"Proteins_plddt.csv"),row.names = F)
pdb_files <-list.files(path=wd1,pattern = ".pdb")
fasta_files=list.files(path=wd2,pattern=".fasta")
pdb_files_no_ext <- tools::file_path_sans_ext(pdb_files)
fasta_files_no_ext <- tools::file_path_sans_ext(fasta_files)

diff_fasta_pdb <- setdiff(fasta_files_no_ext, pdb_files_no_ext)
f_fasta=data.frame(diff_fasta_pdb)
write.csv(f_fasta,str_c(wd1,"missed_proteins.csv"),row.names = F)
missed_fasta=paste0(wd2,diff_fasta_pdb,".fasta")
writeLines(missed_fasta,str_c(wd1,"missed_proteins.txt"))
missed_fasta_1=c("Fasta_Path",missed_fasta)
writeLines(missed_fasta_1,str_c(wd1,"missed_proteins_1.csv"))