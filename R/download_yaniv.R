
#libraries####
library(Biostrings) #not sure if I needed this one
library(rentrez)
library(stringr)

#notes or links####
#https://www.reneshbedre.com/blog/ncbi_sra_toolkit.html
#mention a dot before a funciton for hidden 

#future option
system("esearch -db sra -query 'PRJNA517527' | efetch -format runinfo | grep -v 'Run' | cut -d ',' -f10")
# output
# https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR8504995/SRR8504995.1
# https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR8505002/SRR8505002.1
# https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-15/SRR8504977/SRR8504977.1
# https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR8504976/SRR8504976.1

#data bases options in for future 
entrez_dbs()

# fun <- function() {
#   ANSWER <- readline("Are you a satisfied R user? ")
#   ## a better version would check the answer less cursorily, and
#   ## perhaps re-prompt
#   if (substr(ANSWER, 1, 1) == "n")
#     cat("This is impossible.  YOU LIED!\n")
#   else
#     cat("I knew it.\n")
# }
# if(interactive()) fun()
# }

#functions####
fastq_d_function <- function(runsFile_vector_input){
  
  for(i in 1:length(runsFile_vector_input))
  {
    cat("current input is -",runsFile_vector_input[i])
    
    fastq_d_system <-paste("parallel-fastq-dump --sra-id", runsFile_vector_input[i], "--threads 8 --outdir ./ --split-files --gzip",sep=" ")
    system(fastq_d_system)
  }
  
  # converted into fastq
  system('gunzip *.fastq.gz')
  
}


fun<- function(is_it_example_run){
  
  
  if(is_it_example_run == "example"){
    print("you have selected the example version this will take 5 minutes to complete for id PRJNA605442 ")
    id_input_number <- "PRJNA605442"
  }
  if(is_it_example_run == "notExample"){
    
    #this is a decision you can make -> do you want a manu that provides prefixs or just make it that you must log in id including prefix? 
    # xtrct_idNum_from_PRJNAinputRdLn <- paste(str_extract_all(readline("type in your id number in here (i.e. PRJNA605442) "), "[0-9]", simplify = TRUE), collapse = "")
    # PRJNA_number <- paste("PRJNA",xtrct_idNum_from_PRJNAinputRdLn,sep="")
    
    #this needs the user to type the entire id number include a prefix
    id_input_number <- paste(readline("type in your id number in here (i.e. PRJNA605442 or PRJEB33577 [we should find a small size PRJEB source for example]) "))
    
    #when we make a decision I will need to create try catch to ensure it goes the way we want the user to log in information
    
    
  }
  
  
  id_input_storage <- readline("type in directory name or desired path:\nExample 1: /projectnb2/ct-shbioinf/ykovalski/Workflow/new_name_of_dir\nExample 2: new_name_of_dir")
  
  system(paste("mkdir",id_input_storage,sep = " "))
  
  setwd(id_input_storage)
  
  
  #same thus far now we show an example of a difference between two type of databases. 
  
  
  if (substr(id_input_number, 1, 5) == "PRJEB"){
    cat("This is PRJEB")
    
    input_system_PRJEB_id <- paste("esearch -db bioproject -query",id_input_number,"| efetch -format runinfo -mode xml | xtract -pattern ProjectID -element CenterID", sep = " ") 
    id_insert <- system(input_system_PRJEB_id, intern = TRUE)
    get_ERR_files_names <- paste("esearch -db sra -query",id_insert,"| efetch -format runinfo -mode xml | xtract -pattern SraRunInfo -element Run", sep= " ")
    ERR_vector_preSplit <- system(get_ERR_files_names, intern = TRUE)
    Run_vector<- strsplit(ERR_vector_preSplit, '\t')
  }
  else if(substr(id_input_number, 1, 5) == "PRJNA"){
    
    input_system_PRJNA_id <- paste("esearch -db sra -query",id_input_number,"| efetch -format runinfo | grep -v 'Run' | cut -d ',' -f1> listRunfiles.txt", sep = " ")
    
    system(input_system_PRJNA_id)
    
    Run_reads <- scan("listRunfiles.txt",what="", sep="\n")
    
    
    Run_vector <- as.vector(Run_reads)
    
  }
  
  fastq_d_function(Run_vector)
}


information <- function(){
  #this is a function that holds information elements of the main manu of fun()
  switch(menu(c("general information","example of workflow (it refers to your website but it could be a github link html)","return to main manu"))+1 ,
         cat("done"), paste("help tag in here like the pydamage one --help"), browseURL("https://bio331.devbioinformatics.org/"), main_manu() )
  
}


#only need this if you want options 
main_manu<- function(){
  switch(menu(c("type in your id number with prefix","Run an example with id PRJNA605442","other option"))+1 ,
         cat("done"), fun("notExample"),fun("example"), information() )
}


#after loading the functions and libraries#### 
main_manu()




