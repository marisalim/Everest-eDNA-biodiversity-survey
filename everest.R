#Charles Xu
#March 18 2021
#Everest project with WCS

#load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(hrbrthemes)
library(parallel)
library(MASS)

#set working directory
setwd("C:/Users/Charles/OneDrive - McGill University/mcgill/research/everest")

# unique orders pulled using X reference genomes --------------------------

#read in data
order_genomes_collection<-read.csv("order_genomes_collection_85_clean_nocock.csv")

#number of unique orders pulled using each reference genome
unique_orders<-distinct(order_genomes_collection) %>% group_by(ref) %>% tally

#average number of unique orders pulled using 1-10 reference genomes
for (i in 1:length(unique_orders$ref)){
  sum=0
  com<-combn(unique_orders$ref,i)
  for (j in 1:ncol(com)){
    dist<-distinct(order_genomes_collection) %>%
            filter(ref %in% com[,j]) %>%
            distinct(order) %>%
            count
    sum=sum+dist
  }
  print(sum/ncol(com))
}

#as a check
mean(unique_orders$n) #25.2 is average # of unique orders pulled using 1 genome, which is correct
length(unique(order_genomes_collection$order)) #155 is the total number of unique orders pulled using all 10 genomes, which is correct

#create data frame to plot and add origin
avg_orders<-c(19.66667,36.27778,50.65476,63.42063,75.03175,85.80952,95.97222,105.6667,115)
dat<-data.frame("num_genomes"=seq_along(avg_orders),"avg_orders"=avg_orders)
dat<-rbind(dat,c(0,0))
ggplot(data=dat,aes(x=num_genomes,y=avg_orders)) +
         geom_step() +
         geom_point() +
         scale_x_continuous("number of reference genomes",n.breaks=10,expand=c(0,0)) +
         scale_y_continuous("number of unique orders pulled",n.breaks=25,expand=c(0,0)) +
         ggtitle("Average # of unique orders pulled using X number of reference genomes") +
         geom_line( color="#69b3a2", size=2, alpha=0.9, linetype=2) +
         theme_ipsum() +
         theme(panel.border = element_blank(),panel.grid.minor = element_blank())

#model the asymptotic regression function using non-linear least squares
model<-nls(avg_orders~SSasympOrig(num_genomes,Asym,lrc),data=dat)
summary(model)
pdat<-data.frame(num_genomes=c(0:60))
pdat$fit<-predict(model,newdata=pdat)
plot(avg_orders~num_genomes,data=dat,ylim=c(0,190),xlim=c(0,60))
lines(pdat$num_genomes,pdat$fit,lwd=2)
#add dotted line at ~231, the average maximum number of orders that can be pulled
abline(h=model$m$getPars()["Asym"],lty='dotted')
#add dotted line at ~208, the 90th percentile
abline(h=model$m$getPars()["Asym"]*0.9,lty='dotted')
#add dotted line at ~116, the 50th percentile
abline(h=model$m$getPars()["Asym"]/2,lty='dotted')

#use ggplot to look nice
pdat<-cbind(dat,pdat$fit)
ggplot(data=dat,aes(x=num_genomes,y=avg_orders)) +
  geom_step() +
  geom_point() +
  scale_x_continuous("number of reference genomes",n.breaks=10,expand=c(0,0),limits=c(0,60)) +
  scale_y_continuous("number of unique orders pulled",n.breaks=25,expand=c(0,0),limits=c(0,190)) +
  ggtitle("Average # of unique orders pulled using X number of reference genomes") +
  geom_line(color="#69b3a2",size=2,alpha=0.9,linetype=2) +
  theme_ipsum() +
  theme(panel.border = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15)) +
  geom_line(data=pdat,aes(x=num_genomes,y=fit),size=1,linetype=1) +
  geom_hline(yintercept=model$m$getPars()["Asym"],linetype="dashed",color="red") +
  geom_hline(yintercept=model$m$getPars()["Asym"]*0.9,linetype="dashed",color="red") +
  geom_hline(yintercept=model$m$getPars()["Asym"]/2,linetype="dashed",color="red")


# richness curves ---------------------------------------------------------

#read in data and set order names as row names
filelist<-list.files(pattern='^data\\d+_...\\.csv')
richness_data=list()
for (h in 1:length(filelist)){
  richness_data=c(richness_data,list(data.frame(lapply(filelist,read.csv)[[h]][,-1],
                                   row.names=lapply(filelist,read.csv)[[h]][,1])))
}

#number of unique orders detected per sample
for (m in 1:8){
  print(filelist[m])
  data.frame(richness_data[[m]][1:10]!=0) %>%
    lapply(as.numeric) %>%
    lapply(sum) %>%
    print
}

#This take a long time to run because there are so many unique combinations of 20 samples
#average number of unique orders pulled using 1-10 merged samples
#for each data set
for (j in 1:8){
  avg_list=list()
  #for each number of subset that can be drawn from all samples
  for (i in 1:ncol(richness_data[[j]])){
    total=0
    com_s<-combn(colnames(richness_data[[j]]),i)
    #for each unique combination of samples
    for (k in 1:ncol(com_s)){
      samples=list()
      #combine orders detected
      mclapply(com_s[,k],function(sample){
        samples<<-bind_cols(samples,richness_data[[j]][sample])
      })
      #sum the number of unique orders detected across the combination of samples
      total=total+sum(as.numeric(rowSums(samples)>0))
    }
    #mean number of unique orders detected for a certain subset number
    avg<-mean(total/ncol(combn(colnames(richness_data[[j]]),i)))
    avg_name<-paste("avg",i,sep="")
    avg_list<-rbind(avg_list,cbind(avg_name,avg))
  }
  #for each data set
  data_name<-paste("data",j,sep=".")
  assign(data_name,avg_list)
}

#use ggplot to look nice
lapply(1:8, function(d){
df<-as.data.frame(get(paste("data",d,sep=".")))
df$sample<-1:length(df$avg_name)
df$avg<-as.numeric(unlist(df$avg))
df<-rbind(df,c(0,0,0))
ggplot(data=df,aes(x=sample,y=avg)) +
    geom_step() +
    geom_point() +
    scale_x_continuous("number of samples",breaks=0:50,expand=c(0,0),limits=c(0,50)) +
    scale_y_continuous("number of unique orders",n.breaks=10,expand=c(0,0),limits=c(0,400)) +
    ggtitle("Average # of unique orders detected from X number of samples") +
    geom_line(color="#69b3a2",size=2,alpha=0.9,linetype=2) +
    theme_ipsum() +
    theme(panel.border = element_blank(),panel.grid.minor = element_blank(),
            axis.text.x=element_text(size=15),
            axis.text.y=element_text(size=15)) %>% print
})

#model the asymptotic regression function using non-linear least squares
d=7 #change from 1-8 as appropriately
  df<-as.data.frame(get(paste("data",d,sep=".")))
  df$sample<-1:length(df$avg_name)
  df$avg<-as.numeric(unlist(df$avg))
  df<-rbind(df,c(0,0,0))
  model<-nls(avg~SSasympOrig(sample,Asym,lrc),data=df)
  summary(model)
  pdat<-data.frame(sample=c(0:100))
  pdat$fit<-predict(model,newdata=pdat)
  plot(avg~sample,data=df,ylim=c(0,400),xlim=c(0,54))
  lines(pdat$num_genomes,pdat$fit,lwd=2)
  #add dotted line at ~231, the average maximum number of orders that can be pulled
  abline(h=model$m$getPars()["Asym"],lty='dotted')
  #add dotted line at ~208, the 90th percentile
  abline(h=model$m$getPars()["Asym"]*0.9,lty='dotted')
  #add dotted line at ~116, the 50th percentile
  abline(h=model$m$getPars()["Asym"]/2,lty='dotted')
  
#use ggplot to look nice
  pdat<-cbind(df,pdat$fit)
  ggplot(data=df,aes(x=sample,y=avg)) +
    geom_step() +
    geom_point() +
    scale_x_continuous("number of samples",n.breaks=10,expand=c(0,0),limits=c(0,54)) +
    scale_y_continuous("number of unique orders pulled",n.breaks=25,expand=c(0,0),limits=c(0,150)) +
    ggtitle("Average # of unique orders pulled using X samples") +
    geom_line(color="#69b3a2",size=2,alpha=0.9,linetype=2) +
    theme_ipsum() +
    theme(panel.border = element_blank(),panel.grid.minor = element_blank(),
          axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15)) +
    geom_line(data=pdat,aes(x=sample,y=fit),size=1,linetype=1) +
    geom_hline(yintercept=model$m$getPars()["Asym"],linetype="dashed",color="red") +
    geom_hline(yintercept=model$m$getPars()["Asym"]*0.9,linetype="dashed",color="red") +
    geom_hline(yintercept=model$m$getPars()["Asym"]/2,linetype="dashed",color="red")
  
  