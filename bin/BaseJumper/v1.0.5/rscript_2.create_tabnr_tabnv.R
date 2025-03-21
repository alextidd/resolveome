library(dplyr)
library(magrittr)
library(data.table)

Tab_nv <- data.table::fread("df_nv.tsv",sep="\t",quote = "") %>% 
    reshape2::acast(V2~V1,value.var = "V3",fill = 0)

Tab_nr <- data.table::fread("df_nr.tsv",sep="\t",quote = "") %>% 
    reshape2::acast(V2~V1,value.var = "V3",fill = 0)

pos_nv <- rownames(Tab_nv) %>% strsplit(split = "_") %>% lapply(function(x)paste0(x[1],"_",x[2])) %>% unlist

Tab_nr <- match(pos_nv,rownames(Tab_nr)) %>% Tab_nr[.,]

rownames(Tab_nr) <- rownames(Tab_nv)

Tab_nr <- match(colnames(Tab_nv),colnames(Tab_nr)) %>% Tab_nr[,.]

moutfile_nr <- paste0("Tab_NR.tsv")

moutfile_nv <- paste0("Tab_NV.tsv")

#Write matrices
mheader <- c("",colnames(Tab_nr))
    
write.table(x = mheader %>% t,file = moutfile_nr,append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
  
mheader <- c("",colnames(Tab_nv))

write.table(x = mheader %>% t,file = moutfile_nv,append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

write.table(x = Tab_nr,file = moutfile_nr,append = TRUE,quote = FALSE,sep = "\t",row.names = TRUE,col.names = FALSE)
    
write.table(x = Tab_nv,file = moutfile_nv,append = TRUE,quote = FALSE,sep = "\t",row.names = TRUE,col.names = FALSE)