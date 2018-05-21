rm(list=ls())

#normalize matrix along the row vector
row_nor <- function(input) {
  rows <- dim(input)[1]
  cols <- dim(input)[2]
  output <- matrix(nrow = rows,ncol = cols)
  for (i in 1:rows) {
    s <- sum(input[i,])
    for (j in 1:cols) {
      if (s!=0) {
        output[i,j] <- input[i,j]/s
      } else {
        output[i,j] <- 0
      }
    }
  }
  return(output)
}


#normalize matrix along the column vector
col_nor <- function(input) {
  rows <- dim(input)[1]
  cols <- dim(input)[2]
  output <- matrix(nrow = rows, ncol = cols)
  for (i in 1:rows) {
    s <- sum(input[,i])
    for (j in 1:cols) {
      output[i,j] <- input[i,j]/s
    }
  }
  return(output)
}


#decomposition of atomic relation for odd-length meta-path
split_mat <- function(input) {
  rows <- dim(input)[1]
  cols <- dim(input)[2]
  mid <- sum(input)
  output1 <- matrix(nrow=rows,ncol=mid)
  output1[is.na(output1)] <- 0
  output2 <- matrix(nrow=mid,ncol=cols)
  output2[is.na(output2)] <- 0
  start <- 1
  for (i in 1:rows) {
    num <- sum(input[i,])
    if (num != 0) {
      for (j in start:(start+num-1)) {
        output1[i,j] <- 1
      }
      start <- start + num
    }
  }
  start <- 1
  for (i in 1:cols) {
    num <- sum(input[,i])
    if (num != 0) {
      for (j in start:(start+num-1)) {
        output2[j,i] <- 1
      }
      start <- start + num
    }
  }
  output <- list(output1,output2)
  return(output)
}


#normalization of HeteSim
final_nor <- function(input1, input2) {
  rows1 <- dim(input1)[1]
  cols1 <- dim(input1)[2]
  rows2 <- dim(input2)[1]
  cols2 <- dim(input2)[2]
  output <- matrix(nrow = rows1, ncol = cols2)
  for (i in 1:rows1) {
    for (j in 1:cols2) {
      f1 <- sum(input1[i,]*input2[,j])
      f2 <- norm(as.matrix(input1[i,]),'F')*norm(as.matrix(input2[,j]),'F')
      if (f2!=0) {
        output[i,j] <- f1/f2
      } else {
        output[i,j] <- 0
      }
      
    }
  }
  return(output)
}


#compute HeteSim similarity
hetesim <- function(DG_mat, DD_mat, GG_mat, m_p) {
  drug_num <- dim(DD_mat)[1]
  comp_mat1 <- diag(drug_num)
  comp_mat2 <- diag(drug_num)
  #the length of meta-path is even 
  if (nchar(m_p)%%2==0) {
    #obtain the middle position and then decompose meta-path
    so <- floor(nchar(m_p)/2)-1
    #traverse the first decomposed meta-path
    for (j in 1:so) {
      comp_mat1 <- row_nor(comp_mat1)
      if (substr(m_p,j,j)=='D' && substr(m_p,j+1,j+1)=='D') {
        temp <- row_nor(DD_mat)
        print(j)
      }
      else if (substr(m_p,j,j)=='D' && substr(m_p,j+1,j+1)=='G') {
        temp <- row_nor(DG_mat)
      }
      else if (substr(m_p,j,j)=='G' && substr(m_p,j+1,j+1)=='D') {
        temp <- row_nor(t(DG_mat))
      }
      else if (substr(m_p,j,j)=='G' && substr(m_p,j+1,j+1)=='G') {
        temp <- row_nor(GG_mat)
      }
      comp_mat1 <- comp_mat1 %*% temp
    }
    st <- (nchar(m_p)-1)
    so <- (nchar(m_p)/2+2)-1
    #traverse the second decomposed meta-path
    for (j in st:so) {
      comp_mat2 <- row_nor(comp_mat2)
      if (substr(m_p,j+1,j+1)=='D' && substr(m_p,j,j)=='D') {
        temp <- row_nor(DD_mat)
        print(j)
      }
      else if (substr(m_p,j+1,j+1)=='D' && substr(m_p,j,j)=='G') {
        temp <- row_nor(DG_mat)
      }
      else if (substr(m_p,j+1,j+1)=='G' && substr(m_p,j,j)=='D') {
        temp <- row_nor(t(DG_mat))
      }
      else if (substr(m_p,j+1,j+1)=='G' && substr(m_p,j,j)=='G') {
        temp <- row_nor(GG_mat)
      }
      comp_mat2 <- comp_mat2 %*% temp
    }
    #decomposition of atomic relation
    mid_pos <- nchar(m_p)/2
    if (substr(m_p,mid_pos,mid_pos)=='D' && substr(m_p,mid_pos+1,mid_pos+1)=='D') {
      AB <- split_mat(DD_mat)
    }
    else if (substr(m_p,mid_pos,mid_pos)=='D' && substr(m_p,mid_pos+1,mid_pos+1)=='G') {
      AB <- split_mat(DG_mat)
    }
    else if (substr(m_p,mid_pos,mid_pos)=='G' && substr(m_p,mid_pos+1,mid_pos+1)=='D') {
      AB <- split_mat(t(DG_mat))
    }
    else if (substr(m_p,mid_pos,mid_pos)=='G' && substr(m_p,mid_pos+1,mid_pos+1)=='G') {
      AB <- split_mat(GG_mat)
    }
    A <- AB[[1]]
    B <- AB[[2]]
    comp_mat1 <- row_nor(comp_mat1)
    A <- row_nor(A)
    comp_mat1 <- comp_mat1 %*% A
    comp_mat2 <- row_nor(comp_mat2)
    B <- t(B)
    B <- row_nor(B)
    comp_mat2 <- comp_mat2 %*% B
    #drugsim_mat <- matrix(nrow=drug_num,ncol=drug_num)
    drugsim_mat <- final_nor(comp_mat1,t(comp_mat2))
  } else {
    #the length of meta-path is odd
    #traverse the first decomposed meta-path
    for (j in 1:floor(nchar(m_p)/2)) {
      comp_mat1 <- row_nor(comp_mat1)
      if (substr(m_p,j,j)=='D' && substr(m_p,j+1,j+1)=='D') {
        temp <- row_nor(DD_mat)
        print(1)
      }
      else if (substr(m_p,j,j)=='D' && substr(m_p,j+1,j+1)=='G') {
        temp <- row_nor(DG_mat)
        print(2)
      }
      else if (substr(m_p,j,j)=='G' && substr(m_p,j+1,j+1)=='D') {
        temp <- row_nor(t(DG_mat))
        print(3)
      }
      else if (substr(m_p,j,j)=='G' && substr(m_p,j+1,j+1)=='G') {
        temp <- row_nor(GG_mat)
        print(4)
      }
      comp_mat1 <- comp_mat1 %*% temp
    }
    st <- (nchar(m_p)-1)
    so <- (floor(nchar(m_p)/2+2))-1
    #traverse the second decomposed meta-path
    for (j in st:so) {
      comp_mat2 <- row_nor(comp_mat2)
      if (substr(m_p,j+1,j+1)=='D' && substr(m_p,j,j)=='D') {
        temp <- row_nor(DD_mat)
      }
      else if (substr(m_p,j+1,j+1)=='D' && substr(m_p,j,j)=='G') {
        temp <- row_nor(DG_mat)
      }
      else if (substr(m_p,j+1,j+1)=='G' && substr(m_p,j,j)=='D') {
        temp <- row_nor(t(DG_mat))
      }
      else if (substr(m_p,j+1,j+1)=='G' && substr(m_p,j,j)=='G') {
        temp <- row_nor(GG_mat)
      }
      comp_mat2 <- comp_mat2%*%temp
    }
    #normalization of hetesim
    drugsim_mat <- final_nor(comp_mat1,t(comp_mat2))
  }
  return(drugsim_mat)
}


#corresponds to D^(-1/2), D is a diagonal matrix
fix_no_zero <- function(input) {
  output <- matrix(nrow=dim(input)[1],ncol=dim(input)[2])
  for (i in 1:dim(input)[1]) {
    for (j in 1:dim(input)[2]) {
      if (input[i,j]==0 && i==j) {
        output[i,j] <- 0
      }
      else if (i==j) {
        output[i,j] <- input[i,j]^(-0.5)
      }
      else {
        output[i,j] <- input[i,j]
      }
    }
  }
  return(output)
}


#read data from file
DG_mat <- read.table("E:\\research\\data\\related drug combinations\\drug_final\\DG_mat.txt") #DG_mat represents the relationship between drug and gene
DD_mat_O <- read.table("E:\\research\\data\\related drug combinations\\drug_final\\DD_mat.txt") #DD_mat_O represents the known drug combinations
GG_mat <- read.table("E:\\research\\data\\related drug combinations\\drug_final\\GG_mat.txt") #GG_mat represents the gene-gene interactions
DD_sim_ATC1 <- read.table("E:\\research\\data\\related drug combinations\\drug_final\\DD_sim_ATC1.txt") #DD_sim_ATC1 represents the ATC code similarity at 1-th level
DD_sim_ATC2 <- read.table("E:\\research\\data\\related drug combinations\\drug_final\\DD_sim_ATC2.txt") #DD_sim_ATC2 represents the ATC code similarity at 2-th level
DD_sim_ATC3 <- read.table("E:\\research\\data\\related drug combinations\\drug_final\\DD_sim_ATC3.txt") #DD_sim_ATC3 represents the ATC code similarity at 3-th level
DD_sim_chemical <- read.table("E:\\research\\data\\related drug combinations\\drug_final\\DD_sim_chemical.txt") # DD_sim_chemical represents the chemical structure similarity




DG_mat <- as.matrix(DG_mat)
DD_mat_O <- as.matrix(DD_mat_O)
GG_mat <- as.matrix(GG_mat)
DD_sim_ATC1 <- as.matrix(DD_sim_ATC1)
DD_sim_ATC2 <- as.matrix(DD_sim_ATC2)
DD_sim_ATC3 <- as.matrix(DD_sim_ATC3)
DD_sim_chemical <- as.matrix(DD_sim_chemical)


#obtain drug number in dataset
drug_num <- dim(DD_mat_O)[1]


#calculate drug similarity based on different meta-paths
DD_mat <- DD_mat_O
SDD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DD')
SDDD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DDD')
SDGD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DGD')
SDGDD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DGDD')
SDGGD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DGGD')
SDDDD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DDDD')
SDDGD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DDGD')
SDDDDD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DDDDD')
SDDDGD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DDDGD')
SDDGDD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DDGDD')
SDDGGD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DDGGD')
SDGDDD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DGDDD')
SDGDGD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DGDGD')
SDGGDD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DGGDD')
SDGGGD <- hetesim(DG_mat, DD_mat, GG_mat, m_p='DGGGD')


# attr_num is equaled to that dimension of feature vector plus 1
attr_num <- 20


attr <- data.frame(matrix(nrow=drug_num*(drug_num-1)/2,ncol=attr_num))
attr_pos <- 1


#concatenate feature vector of drug pair
for (m in 2:drug_num) {
  for (n in 1:(m-1)) {
    attr[attr_pos,1] <- DD_sim_ATC1[m,n]
    attr[attr_pos,2] <- DD_sim_ATC2[m,n]
    attr[attr_pos,3] <- DD_sim_ATC3[m,n]
    attr[attr_pos,4] <- DD_sim_chemical[m,n]
    attr[attr_pos,5] <- SDD[m,n] #correspond to meta-path DD
    attr[attr_pos,6] <- 0.5*SDDD[m,n] #correspond to meta-path DDD
    attr[attr_pos,7] <- 0.5*SDGD[m,n] #correspond to meta-path DGD
    attr[attr_pos,8] <- 0.33*SDGDD[m,n] #correspond to meta-path DGDD
    attr[attr_pos,9] <- 0.33*SDGGD[m,n] #correspond to meta-path DGGD
    attr[attr_pos,10] <- 0.33*SDDDD[m,n] #correspond to meta-path DDDD
    attr[attr_pos,11] <- 0.33*SDDGD[m,n] #correspond to meta-path DDGD
    attr[attr_pos,12] <- 0.25*SDDDDD[m,n] #correspond to meta-path DDDDD
    attr[attr_pos,13] <- 0.25*SDDDGD[m,n] #correspond to meta-path DDDGD
    attr[attr_pos,14] <- 0.25*SDDGDD[m,n] #correspond to meta-path DDGDD
    attr[attr_pos,15] <- 0.25*SDDGGD[m,n] #correspond to meta-path DDGGD
    attr[attr_pos,16] <- 0.25*SDGDDD[m,n] #correspond to meta-path DGDDD
    attr[attr_pos,17] <- 0.25*SDGDGD[m,n] #correspond to meta-path DGDGD
    attr[attr_pos,18] <- 0.25*SDGGDD[m,n] #correspond to meta-path DGGDD
    attr[attr_pos,19] <- 0.25*SDGGGD[m,n] #correspond to meta-path DGGGD
    attr[attr_pos,attr_num] <- DD_mat_O[m,n] #the last column denote label
    attr_pos <- attr_pos+1
  }
}

#obtain feature vector of train positive samples trp_set and negative samples n_set 
trp_pos <- 1
n_pos <- 1
trp_set <- data.frame(matrix(0,nrow = 1, ncol = attr_num))
n_set <- data.frame(matrix(0,nrow = 1, ncol = attr_num))
for (j in 1:dim(attr)[1]) {
  if (attr[j,attr_num]==1) {
    trp_set[trp_pos,] <- attr[j,]
    trp_pos <- trp_pos + 1
  }
  else if (attr[j,attr_num]==0) {
    n_set[n_pos,] <- attr[j,]
    n_pos <- n_pos + 1
  } 
}


#the Euclidean distance between negative samples with original point
thr_set_pos <- 1
thr_set <- c()
for (j in 1:dim(n_set)[1]) {
  thr_set[thr_set_pos] <- sum(n_set[j,]^2)
  thr_set_pos <- thr_set_pos + 1
}


#select train negative samples  
trn_pos <- 1
s_t <- sort(thr_set)
threhold <- s_t[dim(trp_set)[1]]
trn_set <- data.frame(matrix(0,1,attr_num))
trn_pos_set <- c()
for (j in 1:dim(n_set)[1]) {
  t_thr <- sum(n_set[j,]^2)
  if (t_thr <= threhold) {
    trn_set[trn_pos,] <- n_set[j,]
    trn_pos_set[trn_pos] <- j
    trn_pos <- trn_pos + 1
  }
}

#construct obtain train set and test set
tr_set <- rbind(trp_set,trn_set)
te_set <- rbind(n_set)


#update order of test samples, in order to combine results from feature-based predictor and network-based predictor  
tep_pos_in_all <- p_in_all[((i-1)*num_te+1):(i*num_te)]
te_set2 <- data.frame(matrix(0,dim(te_set)[1],attr_num))
te_set2_pos <- 1
for (j in 1:dim(attr)[1]) {
  if (attr[j,attr_num]==0) {
    for (m in 1:(attr_num-1)) {
      te_set2[te_set2_pos,m] <- attr[j,m]
      te_set2[te_set2_pos,attr_num] <- 0
    }
    te_set2_pos <- te_set2_pos + 1
  }
}
te_set <- te_set2

feature <- tr_set
test <- te_set
colnames(feature)[attr_num] <- c("label")
colnames(test)[attr_num] <- c("label")


#train SVM  
library(e1071)
sv<-svm(label~.,data=feature,type='C-classification',kernel='radial',probability=TRUE)


#obtain feature-based score based on SVM  
pre<-predict(sv,test,decision.values=TRUE,probability=TRUE)
prob=attr(pre,"probabilities")


# #train NB and obtain feature-based score based on NB
# nb_model <- naiveBayes(label~.,data = feature,laplace=1)
# pred_nb <- predict(nb_model,test,type = "raw")

data_f <- data.frame(prob=prob[,1],obs=test$label)

#network-based base predictor
train_st <- DD_mat
WDD <- train_st
WDG <- DG_mat
WGD <- t(DG_mat)
WGG <- GG_mat


YG <- WGD
YD <- WDD


DD <- rowSums(WDD)
DG <- rowSums(WDG)
GD <- rowSums(WGD)
GG <- rowSums(WGG)

DD <- diag(DD)
DG <- diag(DG)
GD <- diag(GD)
GG <- diag(GG)


DD <- fix_no_zero(DD)
DG <- fix_no_zero(DG)
GG <- fix_no_zero(GG)
GD <- fix_no_zero(GD)

SDD <- DD%*%WDD%*%DD
SDG <- DG%*%WDG%*%GD
SGD <- GD%*%WGD%*%DG
SGG <- GG%*%WGG%*%GG


FDt1 <- YD
FGt1 <- YG
ID <- diag(dim(YD)[1])
IG <- diag(dim(YG)[1])


#the first iteration  
FDt2 <- solve(2*ID+2*(ID-SDD))%*%(SDG%*%FGt1+YD)
FGt2 <- solve(2*IG+2*(IG-SGG))%*%(SGD%*%FDt1+YG)
#it represents threhold
it <- sqrt(sum((FDt2-FDt1)^2))


#iteration
while (it > 0.01) {
  FDt1 <- FDt2
  FGt1 <- FGt2
  FDt2 <- solve(2*ID+2*(ID-SDD))%*%(SDG%*%FGt1+YD)
  FGt2 <- solve(2*IG+2*(IG-SGG))%*%(SGD%*%FDt1+YG)
  it <- sqrt(sum((FDt2-FDt1)^2))
}

FDD <- FDt2

#
probability <- c()
label <- c()
list_pos <- 1
for (m in 2:drug_num) {
  for (n in 1:(m-1)) {
    if (train_st[m,n]==0) {
      probability[list_pos] <- FDD[m,n]
      label[list_pos] <- DD_mat_O[m,n]
      list_pos <- list_pos + 1 
    }
  }
}

data_n <- data.frame(prob=probability,obs=label)


#weighted average ensemble rule(uniform blending)
data_en <- data.frame(f1=data_f$prob, f2=data_n$prob, label=data_n$obs)
a <- sum(data_en$f2)/(sum(data_en$f1)+sum(data_en$f2))
b <- sum(data_en$f1)/(sum(data_en$f1)+sum(data_en$f2))
prob_a <- a*data_en$f1+b*data_en$f2
