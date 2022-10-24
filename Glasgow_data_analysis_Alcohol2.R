
rm(list=ls(all=TRUE))
gc(reset = TRUE)


library(tidyverse)
library(ggplot2)

#-------------------------------------------------------------------
# Function that count the number of friends
#-------------------------------------------------------------------

setwd("C:/Users/Admin/Dropbox/Sophie+Thomas+Stefano/1. Data/GlasgowTeenagers_Data")
setwd("C:/Users/s.nasini/Dropbox/Sophie+Thomas+Stefano/1. Data/GlasgowTeenagers_Data")

# Load and count friendship nominations.

load("Glasgow-friendship.RData")
load("Glasgow-substances.RData")

#Pre-processing on friendships data

# remove  structural absence of the tie ie transform Code 10 into 0
friendship.1[which(friendship.1 == 10)] = 0
friendship.2[which(friendship.2 == 10)] = 0
friendship.3[which(friendship.3 == 10)] = 0

#recoding friendship in increasing degree 
friendship.1[which(friendship.1 == 0)] = 3
friendship.1 = 3 - friendship.1

friendship.2[which(friendship.2 == 0)] = 3
friendship.2 = 3 - friendship.2

friendship.3[which(friendship.3 == 0)] = 3
friendship.3 = 3 - friendship.3


# imputation of NA at each measurement time

# At time t1
NA_t1 = which(is.na(friendship.1),arr.ind = TRUE)

for (i in 1:nrow(NA_t1)) {
  
  if(is.na(friendship.2[NA_t1[i,1], NA_t1[i,2]])){
    
    if(is.na(friendship.3[NA_t1[i,1],NA_t1[i,2]])){ # NA NA NA
      
      friendship.1[NA_t1[i,1],NA_t1[i,2]] = 0
      friendship.2[NA_t1[i,1],NA_t1[i,2]] = 0
      friendship.3[NA_t1[i,1],NA_t1[i,2]] = 0
      
    }else{ # NA NA --
      
      friendship.1[NA_t1[i,1],NA_t1[i,2]] = friendship.3[NA_t1[i,1],NA_t1[i,2]]
      friendship.2[NA_t1[i,1],NA_t1[i,2]] = friendship.3[NA_t1[i,1],NA_t1[i,2]]
      
    } 
    
  }else{ # NA --
    
    friendship.1[NA_t1[i,1],NA_t1[i,2]] = friendship.2[NA_t1[i,1],NA_t1[i,2]]
    
  }
  
}

# At time t2

NA_t2 = which(is.na(friendship.2),arr.ind = T)

for (i in 1:nrow(NA_t2)) {
  
  if(is.na(friendship.3[NA_t2[i,1],NA_t2[i,2]])){ # -- NA NA
    
    friendship.2[NA_t2[i,1],NA_t2[i,2]] = friendship.1[NA_t2[i,1],NA_t2[i,2]]
    friendship.3[NA_t2[i,1],NA_t2[i,2]] = friendship.1[NA_t2[i,1],NA_t2[i,2]]
    
  }else{ # -- NA --
    
    friendship.2[NA_t2[i,1],NA_t2[i,2]] = round(0.5*(friendship.1[NA_t2[i,1],NA_t2[i,2]] + friendship.3[NA_t2[i,1],NA_t2[i,2]]))
    
  }
  
}

# At time t3

NA_t3 = which(is.na(friendship.3),arr.ind = T)

for (i in 1:nrow(NA_t3)) {
  
  friendship.3[NA_t3[i,1],NA_t3[i,2]] = friendship.2[NA_t3[i,1],NA_t3[i,2]]
  
}


# Counting the degree of friendship of each child at each measurement time

d_t1 = as.vector(rowSums(friendship.1))
d_t2 = as.vector(rowSums(friendship.2))
d_t3 = as.vector(rowSums(friendship.3))


#Imputation of NA in consumption at each measurement time

alcohol = as.data.frame(alcohol)

IND_d_t1 = which(d_t1 == 0)
IND_d_t2 = which(d_t2 == 0)
IND_d_t3 = which(d_t3 == 0)

AM_1 = friendship.1[-IND_d_t1,-IND_d_t1]
AM_2 = friendship.2[-IND_d_t2,-IND_d_t2]
AM_3 = friendship.3[-IND_d_t3,-IND_d_t3]

d_t1 = d_t1[-IND_d_t1]
d_t2 = d_t2[-IND_d_t2]
d_t3 = d_t3[-IND_d_t3]

X01 = alcohol[,1]
x01 = alcohol[-IND_d_t1,1]

X02 = alcohol[,2]
x02 = alcohol[-IND_d_t2,2]

X03 = alcohol[,3]
x03 = alcohol[-IND_d_t3,3]



#-------------------------------------------------------------------
# Inputation at time t1
#-------------------------------------------------------------------

IND1_NA = which(is.na(x01))
IND1_NNA = which(!is.na(x01))
X1_NNA = x01[IND1_NNA]

if(length(IND1_NA)!=0){
  
  x01[IND1_NA] = mean(x01, na.rm = T)
  for (i in 1:5) {
    
    x01 = round((AM_1%*%x01)/d_t1)
    
    for(j in IND1_NA){
      if(!is.na(x02[j])){
        x01[j] = round(0.5*(x01[j] + X02[j]))
      }else{
        if(!is.na(X03[j])){
          x01[j] = round(0.5*(x01[j] + X03[j]))
        }
      }
    }
    
    x01[IND1_NNA] = X1_NNA
    
  }
  
  #cbind(X01[-IND_d_t1], x01)
  
  X01[-IND_d_t1] = x01
  
  X01[which(is.na(X01))] = round(mean(X01, na.rm = T))
  
  for(j in which(is.na(X01))){
    if(!is.na(X02[j])){
      X01[j] = round(0.5*(X01[j] + X02[j]))
    }else{
      if(!is.na(X03[j])){
        X01[j] = round(0.5*(X01[j] + X03[j]))
      }
    }
  }
}  

#cbind(alcohol[,1],X01)
alcohol[,1] = X01

#-------------------------------------------------------------------
# Inputation at time t2
#-------------------------------------------------------------------

IND2_NA = which(is.na(x02))
IND2_NNA = which(!is.na(x02))
X2_NNA = x02[IND2_NNA]

if(length(IND2_NA)!=0){
  x02[IND2_NA] = round(mean(x02, na.rm = T))
  
  for (i in 1:5) {
    
    x02 = round((AM_2%*%x02)/d_t2)
    
    for(j in IND2_NA){
      if(!is.na(x03[j])){
        x02[j] = round(0.5*(x02[j] + x03[j]))
      }
    }
    
    x02[IND2_NNA] = X2_NNA
    
  }
  
  #cbind(X02[-IND_d_t2], x02)
  
  X02[-IND_d_t2] = x02
  
  X02[which(is.na(X02))] = round(mean(X02, na.rm = T))
  
  for(j in which(is.na(X02))){
    if(!is.na(X03[j])){
      X02[j] = 0.5*(X02[j] + X03[j])
    }
  }
  
}

#cbind(alcohol[,2],X02)

alcohol[,2] = X02

#-------------------------------------------------------------------
# Inputation at time t3
#-------------------------------------------------------------------

IND3_NA = which(is.na(x03))
IND3_NNA = which(!is.na(x03))
X3_NNA = x03[IND3_NNA]

if(length(IND3_NA)!=0){
  
  x03[IND3_NA] = round(mean(x03, na.rm = T))
  
  for (i in 1:5) {
    
    x03 = round((AM_3%*%x03)/d_t3)
    
    x03[IND3_NNA] = X3_NNA
    
  }
  
  #cbind(X03[-IND_d_t3], x03)
  
  X03[-IND_d_t3] = x03
  
  X03[which(is.na(X03))] =round( mean(X03, na.rm = T))
}

#cbind(alcohol[,3],X03)
alcohol[,3] = X03

# rescaling consumption degree from 0 to 3
alcohol = alcohol - 1

#-------------------------------------------------------------------------
# Matrix of differents combinations of three elements of a set
#-------------------------------------------------------------------------

uplet3 = function(set){
  
  M = NULL
  for (i in set) {
    for (j in set) {
      for (k in set) {
        M = rbind(M,c(i,j,k))
      }
      
    }
  } 
  return(M) 
  
}



#------------------------------------------
# data transformation
#-------------------------------------------

#counting the occurence of each combination in the dataset of alcohol

m = uplet3(0:4)
ocurence = rep(0,nrow(m))

for (i  in 1:nrow(m)) {
  for (j in 1:nrow(alcohol)) {
    if(sum(m[i,]==alcohol[j,])==3){ocurence[i] = ocurence[i]+1}
  }
}

alcohol2 = cbind.data.frame(m,ocurence)
alcohol2 = alcohol2[alcohol2[,4]!=0,]


#plots of frenquencies of each entries

setwd("C:/Users/Admin/Dropbox/Sophie+Thomas+Stefano/5. Results/PLOTS")
setwd("C:/Users/s.nasini/Dropbox/Sophie+Thomas+Stefano/5. Results/PLOTS")

filename = paste("BarPlotAdoptionAlcohol.pdf")
pdf(filename )

Labels = rep('a', nrow(alcohol2))

OR = order(alcohol2[,4], decreasing = TRUE)

for(i in 1:nrow(alcohol2)){
  
  Labels[i] = paste(alcohol2[OR[i], 1], alcohol2[OR[i], 2],alcohol2[OR[i], 3])
}
barplot(sort(alcohol2[,4], decreasing=T), axes = T, space = 2, horiz = TRUE, xlim = c(0,100), names.arg = Labels, cex.names = 1.1, xlab= "frequencies",las=1, cex = 1.1, main = "")
dev.off()

# Transformation in binary entries 

#dataframe with same size as alcohol2
alcohol2.b1 = alcohol2
alcohol2.b1[alcohol2.b1==1] = 1
alcohol2.b1[alcohol2.b1==2] = 1
alcohol2.b1[alcohol2.b1==3] = 1
alcohol2.b1[alcohol2.b1==4] = 1
alcohol2.b1[,4] = alcohol2[,4]

# reduced dataframe 
m.b = uplet3(0:1)
ocurence = rep(0,nrow(m.b))

for (i  in 1:nrow(m.b)) {
  for (j in 1:nrow(alcohol2.b1)) {
    if(sum(m.b[i,]==alcohol2.b1[j,])==3){ocurence[i] = ocurence[i]+alcohol2.b1[j,4]}
  }
}

alcohol2.b2 = cbind.data.frame(m.b,ocurence)
alcohol2.b2 = alcohol2.b2[alcohol2.b2[,4]!=0,]

#  consistent entries with the model when 1 is not an adopter
index1 = NULL

for (i in 1:nrow(alcohol2.b2)) {
  if( ( (alcohol2.b2[i, 1] <= alcohol2.b2[i, 2]) & (alcohol2.b2[i, 2]<=alcohol2.b2[i, 3]))){
    index1 = c(index1,i)
  }
}

alcohol2.b2 = alcohol2.b2[index1, ]

# Plot frequencies of new coded data

setwd("C:/Users/Admin/Dropbox/Sophie+Thomas+Stefano/5. Results/PLOTS")
setwd("C:/Users/s.nasini/Dropbox/Sophie+Thomas+Stefano/5. Results/PLOTS")

filename = paste("BarPlotBinaryAdoptionAlcohol.pdf")
pdf(filename )
Labels = rep('a', nrow(alcohol2.b2))

OR = order(alcohol2.b2[,4], decreasing = TRUE)

for(i in 1:nrow(alcohol2.b2)){
  
  Labels[i] = paste(alcohol2.b2[OR[i], 1], alcohol2.b2[OR[i], 2],alcohol2.b2[OR[i], 3])
  
}
barplot(sort(alcohol2.b2[,4], decreasing=T), axes = T, space = 0.5, horiz=TRUE,xlim = c(0,65),names.arg = Labels, cex.names = 1.1, cex = 1.1, xlab= "frequencies",las=1)

dev.off()


#Counting entries consistent with the model when 1 is not an adopter
index = NULL
freq = 0

for (i in 1:nrow(alcohol2.b1)) {
  if(((alcohol2.b1[i, 1] <= alcohol2.b1[i, 2]) & (alcohol2.b1[i, 2]<=alcohol2.b1[i, 3]))){
    index = c(index,i)
    freq = freq + alcohol2.b1[i, 4]
  }
}

#data set of inconsistent data
uncons_alc = alcohol2[-index,]

#indexes of inconsistent data in alcohol
Inc_Ind = NULL
for (i  in 1:nrow(uncons_alc)) {
  for (j in 1:nrow(alcohol)) {
    if(sum(uncons_alc[i,1:3]==alcohol[j,])==3){Inc_Ind = c(Inc_Ind,j)} 
  }
}

#-----------------------------------------------------------------------------------------------------------
# Splitting data into influencials and imitators
#------------------------------------------------------------------------------------------------------------
d_t1 = as.vector(rowSums(friendship.1))
d_t2 = as.vector(rowSums(friendship.2))
d_t3 = as.vector(rowSums(friendship.3))

H = round((d_t1+d_t2+d_t3)/3)
#a = c(0.5, 0.6,  0.7, 0.8, 0.9)

a = c(0.5,  0.7,  0.9)


#-------------------------------------------------------------------------------
H2 = H[-Inc_Ind]
Partition1 = which(H2 > quantile(H2,probs = a[1]))
Partition2 = which(H2 > quantile(H2,probs = a[2]))
Partition3 = which(H2 > quantile(H2,probs = a[3]))
#Partition4 = which(H2 >= quantile(H2,probs = a[4]))
#Partition5 = which(H2 >= quantile(H2,probs = a[5]))

#---------------------------------------------------------------------------------
alcohol1 = alcohol[-Inc_Ind,]
alcohol1[alcohol1==1] = 1
alcohol1[alcohol1==2] = 1
alcohol1[alcohol1==3] = 1
alcohol1[alcohol1==4] = 1

#Matrix of adoptions with repect to partitions

x_par12 = matrix(0,3,2)
#Number of influentials and imitators in the partition
M_par12 = c(length(H2[Partition1]),length(H2[-Partition1]))

#population at time t1
pop_inf_par1_t1 = alcohol1$t1[Partition1]
pop_imi_par1_t1 = alcohol1$t1[-Partition1]
# adoption at time t1
x_par12[1,1] = sum(pop_inf_par1_t1)
x_par12[1,2] = sum(pop_imi_par1_t1)
# population of non adopters at time t1
pop_inf_par1_t2 =alcohol1$t2[Partition1][pop_inf_par1_t1==0]
pop_imi_par1_t2 = alcohol1$t2[-Partition1][pop_imi_par1_t1==0]
#adoptions at time t2
x_par12[2,1] = sum(pop_inf_par1_t2)
x_par12[2,2] = sum(pop_imi_par1_t2)
# population of non adopters at time t2
pop_inf_par1_t3 =alcohol1$t3[Partition1][pop_inf_par1_t1==0 & alcohol1$t2[Partition1]==0]
pop_imi_par1_t3 = alcohol1$t3[-Partition1][pop_imi_par1_t1==0 & alcohol1$t2[-Partition1]==0]
#adoptions at time t3
x_par12[3,1] = sum(pop_inf_par1_t3)
x_par12[3,2] = sum(pop_imi_par1_t3)


#---------------------------------------------------------------------------------
x_par22 = matrix(0,3,2)


#Number of influentials and imitators in the partition
M_par22 = c(length(H2[Partition2]),length(H2[-Partition2]))

#population at time t1
pop_inf_par2_t1 = alcohol1$t1[Partition2]
pop_imi_par2_t1 = alcohol1$t1[-Partition2]
# adoption at time t1
x_par22[1,1] = sum(pop_inf_par2_t1)
x_par22[1,2] = sum(pop_imi_par2_t1)
# population of non adopters at time t1
pop_inf_par2_t2 =alcohol1$t2[Partition2][pop_inf_par2_t1==0]
pop_imi_par2_t2 = alcohol1$t2[-Partition2][pop_imi_par2_t1==0]
#adoptions at time t2
x_par22[2,1] = sum(pop_inf_par2_t2)
x_par22[2,2] = sum(pop_imi_par2_t2)
# population of non adopters at time t2
pop_inf_par2_t3 =alcohol1$t3[Partition2][pop_inf_par2_t1==0 & alcohol1$t2[Partition2]==0]
pop_imi_par2_t3 = alcohol1$t3[-Partition2][pop_imi_par2_t1==0 & alcohol1$t2[-Partition2]==0]
#adoptions at time t3
x_par22[3,1] = sum(pop_inf_par2_t3)
x_par22[3,2] = sum(pop_imi_par2_t3)

#---------------------------------------------------------------------------------
x_par32 = matrix(0,3,2)
#Number of influentials and imitators in the partition
M_par32 = c(length(H2[Partition3]),length(H2[-Partition3]))

#population at time t1
pop_inf_par3_t1 = alcohol1$t1[Partition3]
pop_imi_par3_t1 = alcohol1$t1[-Partition3]
# adoption at time t1
x_par32[1,1] = sum(pop_inf_par3_t1)
x_par32[1,2] = sum(pop_imi_par3_t1)
# population of non adopters at time t1
pop_inf_par3_t2 =alcohol1$t2[Partition3][pop_inf_par3_t1==0]
pop_imi_par3_t2 = alcohol1$t2[-Partition3][pop_imi_par3_t1==0]
#adoptions at time t2
x_par32[2,1] = sum(pop_inf_par3_t2)
x_par32[2,2] = sum(pop_imi_par3_t2)
# population of non adopters at time t2
pop_inf_par3_t3 =alcohol1$t3[Partition3][pop_inf_par3_t1==0 & alcohol1$t2[Partition3]==0]
pop_imi_par3_t3 = alcohol1$t3[-Partition3][pop_imi_par3_t1==0 & alcohol1$t2[-Partition3]==0]
#adoptions at time t3
x_par32[3,1] = sum(pop_inf_par3_t3)
x_par32[3,2] = sum(pop_imi_par3_t3)

setwd("C:/Users/Admin/Dropbox/Sophie+Thomas+Stefano/1. Data")

save(x_par12, x_par22, x_par32, M_par12, M_par22, M_par32,
     file = "alcohol_data_ready_for_estimation2.RData")



