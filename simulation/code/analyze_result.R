
################# making plot for fixed effect, random effect covariance, correlation matrix as a whole ####
rm(list=ls())
library("tikzDevice")
library("ggplot2")

setwd("~/simulation/result")
# each res stores one file: {dls, newdls, XL_dls}*{ind, dep250, dep500}*{20,50,100,150,200}
res <- vector("list",45)


# read file into res 
method <- c("dls","newdls","XL_dls")
subset <- c(20,50,100,150,200)
partition <- c("ind", "dep250", "dep500")
l <- list(m = 1:3, s = 1:5, p = 1:3)
comb <- expand.grid(l)
for(i in 1:45){
  
  m <- as.numeric(comb[i,][1])
  s <- as.numeric(comb[i,][2])
  p <- as.numeric(comb[i,][3])
  data_name <- paste0("sim_",method[m],"_k_",subset[s],"_accuracy_",partition[p],".rds")
  res[[i]] <- readRDS(data_name)[[1]]
  
}

# in total 12 variables

accur_para <- c(  "accur_beta_2D","accur_beta_3D","accur_beta_4D",
                  "accur_dmat_2D","accur_dmat_3D","accur_dmat_4D","accur_dmat_5D","accur_dmat_6D",
                  "accur_Rmat_2D","accur_Rmat_3D")

# each stores one accur_para 
accur_list <- vector("list", 10)
comb$m <- ifelse(comb$m==1, "LS-WASP", ifelse(comb$m==2,"LS-WASP-NEW","DPMC"))
comb$s <- ifelse(comb$s==1, "20", ifelse(comb$s==2, "50", ifelse(comb$s==3, "100",  ifelse(comb$s==4, "150", "200"))))
comb$p <- ifelse(comb$p==1, "ind", ifelse(comb$p==2, "dep250",  "dep500"))


comb_extend <- comb[rep(seq_len(nrow(comb)), each=10),]

for(i in 1:10){
  
  # stock up the twelve files for each var
  column_value <-  unlist(lapply(res, function(x) x$accuracy_matrix[,accur_para[i]]))
  column_value <- matrix(column_value, ncol=1)
  
  accur_list[[i]] <- cbind(column_value, comb_extend)
}


accur_beta_2D <- accur_list[[1]]
accur_beta_3D <- accur_list[[2]]
accur_beta_4D <- accur_list[[3]] 

accur_dmat_2D <- accur_list[[4]]
accur_dmat_3D <- accur_list[[5]] 
accur_dmat_4D <- accur_list[[6]] 
accur_dmat_5D <- accur_list[[7]] 
accur_dmat_6D <- accur_list[[8]] 


accur_Rmat_2D_disjoint <- accur_list[[9]]
accur_Rmat_3D_disjoint <- accur_list[[10]]

accur_Rmat_2D <- accur_list[[9]]
accur_Rmat_3D <- accur_list[[10]]



#################### Make the tables #########################

accur_list <- lapply(accur_list, function(x) x=x[x$m!="LS-WASP",])
accur_list <- lapply(accur_list, function(x) {x$m[x$m=="LS-WASP-NEW"]="LS-WASP"
return(x)})


accur_beta_2D <- accur_list[[1]]
accur_beta_3D <- accur_list[[2]]
accur_beta_4D <- accur_list[[3]] 

accur_dmat_2D <- accur_list[[4]]
accur_dmat_3D <- accur_list[[5]] 
accur_dmat_4D <- accur_list[[6]] 
accur_dmat_5D <- accur_list[[7]] 
accur_dmat_6D <- accur_list[[8]] 


accur_Rmat_2D_disjoint <- accur_list[[9]]
accur_Rmat_3D_disjoint <- accur_list[[10]]

accur_Rmat_2D <- accur_list[[9]]
accur_Rmat_3D <- accur_list[[10]]

accur_list <- list(accur_beta_2D,accur_beta_3D,accur_beta_4D,
                   accur_dmat_2D,accur_dmat_3D,accur_dmat_4D,accur_dmat_5D,accur_dmat_6D,
                   accur_Rmat_2D_disjoint,accur_Rmat_3D_disjoint,accur_Rmat_2D,accur_Rmat_3D)
summary_list <- vector("list",length(accur_list))
names(summary_list) <- c("sum_beta_2D","sum_beta_3D","sum_beta_4D",
              "sum_dmat_2D","sum_dmat_3D","sum_dmat_4D","sum_dmat_5D","sum_dmat_6D",
              "sum_Rmat_2D_disjoint","sum_Rmat_3D_disjoint","sum_Rmat_2D","sum_Rmat_3D")
sum_ind <- 1
for(accur_temp in accur_list){
  
  summary_temp <- aggregate(accur_temp[accur_temp$s!=20,"column_value"], by=list(accur_temp[accur_temp$s!=20,"m"],accur_temp[accur_temp$s!=20,"s"],accur_temp[accur_temp$s!=20,"p"]),
                            function(x) c(mean(x), sd(x)))
  summary_temp <- cbind(summary_temp[,1:3],summary_temp$x)
  colnames(summary_temp) <- c("m","s","p","mean","sd")
  summary_temp <- reshape(summary_temp, idvar = c("m","p"), timevar ="s", direction = "wide")
  summary_temp <- summary_temp[,c("m","p","mean.50","sd.50","mean.100","sd.100","mean.150","sd.150","mean.200","sd.200")]
  summary_list[[sum_ind]] <- summary_temp
  sum_ind <- sum_ind + 1 
}

library("xtable")
## make the table 
sum_table <- function(sum_name_list, sum_dat_list, sum_ind_overlap, partition){
  
  k <- 1
  mean_list <- c("m",sapply(partition, function(x) paste0("mean.",x)))
  sd_list <- c("m",sapply(partition, function(x) paste0("sd.",x)))
  summary_list <- sum_dat_list 
  
  for(i in sum_name_list){
  
  sum_mean <- summary_list[[i]][summary_list[[i]]$p==sum_ind_overlap,mean_list]
  sum_se <- summary_list[[i]][summary_list[[i]]$p==sum_ind_overlap,sd_list]
  sum_mean <- as.matrix(sum_mean)
  sum_mean <- cbind(sum_mean[,1], apply(sum_mean[,2:ncol(sum_mean)],2, function(x) round(as.numeric(x),3) ))
  sum_se <- as.matrix(sum_se)
  sum_se <- cbind(sum_se[1], apply(sum_se[,2:ncol(sum_mean)],2, function(x) round(as.numeric(x),3) ))
  
  M <- matrix(as.vector(rbind(as.character(sum_mean),
                              paste("(",as.vector(sum_se),")", sep=""))), nrow=4)
  M[c(2,4),1] <- c("","")
  colnames(M) <- c("Method",partition)
  
  if(k!=1){
    
    M <- M[,-1]
  }
 
  print(xtable(M),type="latex",include.rownames=FALSE)
  k <- k + 1
}
}



sum_list <- names(summary_list)
## according to the order in the paper 
# beta - ind 
sum_name_list <- sum_list[1:3]
sum_dat_list <- summary_list
sum_ind_overlap <- "ind"
partition <- c(50,100,150,200)

sum_table(sum_name_list, sum_dat_list, sum_ind_overlap, partition)

# D - ind
sum_name_list <- sum_list[4:8]
sum_dat_list <- summary_list
sum_ind_overlap <- "ind"
partition <- c(50,100,150,200)

sum_table(sum_name_list, sum_dat_list, sum_ind_overlap, partition)


# beta - dep
sum_name_list <- sum_list[1:3]
sum_dat_list <- summary_list
sum_ind_overlap <- "dep250"
partition <- c(50,100,150,200)

sum_table(sum_name_list, sum_dat_list, sum_ind_overlap, partition)


sum_name_list <- sum_list[1:3]
sum_dat_list <- summary_list
sum_ind_overlap <- "dep500"
partition <- c(50,100,150,200)

sum_table(sum_name_list, sum_dat_list, sum_ind_overlap, partition)


# D - dep
sum_name_list <- sum_list[4:8]
sum_dat_list <- summary_list
sum_ind_overlap <- "dep250"
partition <- c(50,100,150,200)

sum_table(sum_name_list, sum_dat_list, sum_ind_overlap, partition)


sum_name_list <- sum_list[4:8]
sum_dat_list <- summary_list
sum_ind_overlap <- "dep500"
partition <- c(50,100,150,200)

sum_table(sum_name_list, sum_dat_list, sum_ind_overlap, partition)

# R - dep
sum_name_list <- sum_list[9:10]
sum_dat_list <- summary_list
sum_ind_overlap <- "ind"
partition <- c(50,100,150,200)

sum_table(sum_name_list, sum_dat_list, sum_ind_overlap, partition)


sum_name_list <- sum_list[11:12]
sum_dat_list <- summary_list
sum_ind_overlap <- "dep250"
partition <- c(50,100,150,200)

sum_table(sum_name_list, sum_dat_list, sum_ind_overlap, partition)


sum_name_list <- sum_list[11:12]
sum_dat_list <- summary_list
sum_ind_overlap <- "dep500"
partition <- c(50,100,150,200)

sum_table(sum_name_list, sum_dat_list, sum_ind_overlap, partition)

# to produce the table we want we have to manually 
# 1. only need the part from \begin{tabular} to \end{tabular}
# 1. change \begin{tabular}{llllll} to p{1.35cm}p{0.7cm}p{0.7cm}p{0.7cm} p{0.7cm}p{0.7cm} for the first table
#                                   to p{0.7cm}p{0.7cm}p{0.7cm} p{0.7cm}p{0.7cm} for the rest 
# 2. \begin{table}[ht]
   #%\centering
    #\footnotesize
    #\begin{spacing}{1.2}\centering
    #\begin{tabular}{c@{\quad}c@{\quad}c}
    #$2D$   &   $3D$  &   $4D$ \\

    #\end{tabular}
    #\end{spacing}
    #\end{table}


#p_list <- vector("list",length(accur_list))
#names(p_list) <- c("p_beta_2D","p_beta_3D","p_beta_4D",
#                         "p_dmat_2D","p_dmat_3D","p_dmat_4D","p_dmat_5D","p_dmat_6D",
#                         "p_Rmat_2D_disjoint","p_Rmat_3D_disjoint","p_Rmat_2D","p_Rmat_3D")
#p_ind <- 1
#for(accur_temp in accur_list){

#  p_temp <- aggregate(accur_temp[,"column_value"], by=list(accur_temp[,"s"],accur_temp[,"p"]), 
#                                            function(x) t.test(x[accur_temp$m=="LS-WASP"], x[accur_temp$m=="DPMC"])$p.value)
#  colnames(p_temp) <- c("s","p","p_value")
#  p_list[[p_ind]] <- p_temp
#  p_ind <- p_ind + 1 
#}

