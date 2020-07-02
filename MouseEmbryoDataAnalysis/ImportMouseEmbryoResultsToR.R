#Code to read in the alevin quantification data for the mouse embryo data, both using EM and not using EM
  #Note that to use tximeta you will have to copy the aux_info folder from the EM results alevin folder to the NoEM results folder
  #This information will not be used but tximeta is expecting it to be there to work

#Specify libPaths to ensure the correct one is used
.libPaths("/nas/longleaf/home/skvanbur/lib64/R/library")
UseTxiMeta <- TRUE
if(UseTxiMeta==TRUE){
  library(tximeta)
}else{
  library(tximport)
}
library(fishpond)

base_dirT <- "/pine/scr/s/k/skvanbur/SingleCellProject/MouseEmbryoData/"

for(j in 1:2){
  if(j==1){
    UseEM <- TRUE
  }else if(j==2){
    UseEM <- FALSE
  }
  
  if(UseEM==TRUE){
    base_dir <- file.path(base_dirT, "UsingEM/")
    fil_mod <- ""
    save_mod <- "UsingEM"
    
  }else{
    base_dir <- file.path(base_dirT, "NoEM/")
    fil_mod <- "_noem"
    save_mod <- "NoEM"
  }
  
  if(UseTxiMeta==TRUE){
    save_mod2 <- "TximetaObj"
  }else{
    save_mod2 <- ""
  }
  
  SaveDir <- file.path(base_dirT)
  if(!dir.exists(SaveDir)){dir.create(SaveDir, recursive = T)}
  for(val in 1:3){
    if(val==1){
      num <- "8.0"
    }else if(val==2){
      num <- "8.25"
    }else if(val==3){
      num <- "8.5"
    }
    
    curr_fil <- file.path(base_dir, paste0(num, fil_mod), "alevin", "quants_mat.gz")
    
    
    File_To_Import <- data.frame(curr_fil, stringsAsFactors = F)
    colnames(File_To_Import) <- "files"
    
    File_To_Import$names <- as.character(val)
    if(UseTxiMeta==TRUE){
      #setTximetaBFC(base_dirT)
      assign(paste0("QuantObj", num, save_mod), tximeta(File_To_Import, type = "alevin", alevinArgs = list(filterBarcodes = FALSE, tierImport = TRUE)))
    }else{
      assign(paste0("QuantObj", num, save_mod), tximport(File_To_Import[1,1], type = "alevin", alevinArgs = list(filterBarcodes = FALSE, tierImport = TRUE)))
    }


    
    print(paste0("Current file that was just imported is ", File_To_Import[1,1]))
    print(paste0("Object Name is ", paste0("QuantObj", num, save_mod)))
    print(paste0("Current save file is ", paste0(SaveDir, "QuantObj", num, save_mod, save_mod2,  ".RData")))
    save(list = c(paste0("QuantObj", num, save_mod)), file = paste0(SaveDir, "QuantObj", num, save_mod, save_mod2, ".RData"))
    
    rm(list = paste0("QuantObj", num, save_mod))
    print(gc())
  }
}


