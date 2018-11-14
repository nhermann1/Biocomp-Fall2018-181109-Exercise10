setwd("~/Documents/Biocomp/Biocomp-Fall2018-181109-Exercise10/")
rm(list = ls())

#Question 1
#To test nested models, treat them with a null and alternate hypothesis
#Comparing with null is detailed in lecture notes
#For this, the null is that they're equal and that you should take the simpler, linear model
#Be sure to get data.txt
Data<-read.table("data.txt",header = TRUE,sep = ",") #get data file to test models on

linearNLL<-function(p,x,y){
  a=p[1] 
  b=p[2] 
  sigma=exp(p[3])
  expected=a+x*b
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll) 
} #custom function with a linear expected model

quadNLL<-function(p,x,y){
  a=p[1]
  b=p[2]
  c=p[3]
  sigma=exp(p[4])
  expected=a+b*x+c*x^2
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll) 
} #custom function with a quadratic expected model

initialGuess=c(1,1,1)
linearFIT=optim(par=initialGuess,fn=linearNLL,x=Data$x,y=Data$y)
print(linearFIT) #Gives all the parms and likelihood for linear model when optimized

initialGuess=c(1,1,1,1) 
quadFIT=optim(par=initialGuess,fn=quadNLL,x=Data$x,y=Data$y)
print(quadFIT) #Gives all the parms and likelihood for quadratic model when optimized

test=2*(linearFIT$value-quadFIT$value) #test statistic for likelihood
df=length(quadFIT$par)-length(linearFIT$par) #degrees of freedom for these models

1-pchisq(test,df) #Do 1 minus of chi-squared p-value, so output is likelihood of accepting alternate hypothesis
#Essentially a zero-chance to accept the null that they're equally likely models
#Looking at the values, suggests that the quadratic model is more likely

#Question 2
#I think it wants output from a three scenarios of alpha values
#Most efficient if one is the first logical statement is false
  #Then second is false
  #Then both are true, which should lead to coexistence
#Model should predict the population of both species at various time steps
library(deSolve) #Need for the ode function to use the custom function
library(ggplot2) #Need for good plots
diffSolver<-function(t,y,p) {
  N1=y[1]
  N2=y[2]
  R1=p[1]
  K1=p[2]
  alpha12=p[3]
  R2=p[4]
  K2=p[5]
  alpha21=p[6]
  dN1dt=R1*(1-(N1+alpha12*N2))*N1 #As written, alpha11 is 1 always
  dN2dt=R2*(1-(N2+alpha21*N1))*N2 #As written, alpha22 is 1 always
  return(list(c(dN1dt,dN2dt)))
} #custom function with the two equations to solve 

#Case 1
p2champ=c(0.2,100,1.5,0.2,100,0.5) #ordered R,K,alpha for 1 then for 2--N1 will be extirpated because the competitive effect of 2 is greater than of 1 on itself (i.e. 1)
y=c(0.5,0.5) #Starting points for N1 and N2
t=c(1:100) #Number of "time steps" to take

sim2champ=ode(y=y,times=t,func=diffSolver,parms=p2champ) #run iterations to solve these equations
Output2champ<-data.frame(time=sim2champ[,1],N1=sim2champ[,2],N2=sim2champ[,3]) #insert to data frame

ggplot(Output2champ,aes(x=time))+
  geom_point(aes(y=N1),color="red")+
  geom_point(aes(y=N2),color="green")+
  ylim(0,1)

#Case 2
p1champ=c(0.25,100,0.5,0.25,100,1.5) #ordered R,K,alpha for 1 then for 2--2 is extirpated with lower competitive effect

sim1champ=ode(y=y,times=t,func=diffSolver,parms=p1champ)
Output1champ<-data.frame(time=sim1champ[,1],N1=sim1champ[,2],N2=sim1champ[,3])

ggplot(Output1champ,aes(x=time))+
  geom_point(aes(y=N1),color="red")+
  geom_point(aes(y=N2),color="green")+
  ylim(0,1)

#Case 3
pEQcoexist=c(0.15,100,0.5,0.15,100,0.5) #ordered R,K,alpha for 1 then for 2--both are equal so they will go to equal pop sizes 

simEQcoexist=ode(y=y,times=t,func=diffSolver,parms=pEQcoexist)
OutputEQcoexist<-data.frame(time=simEQcoexist[,1],N1=simEQcoexist[,2],N2=simEQcoexist[,3])

ggplot(OutputEQcoexist,aes(x=time))+
  geom_point(aes(y=N1),color="red")+
  geom_point(aes(y=N2),color="green")+
  ylim(0,1)

#Case 4
pDIFFcoexist=c(0.2,100,0.5,0.2,100,0.25) #ordered R,K,alpha for 1 then for 2--whichever has the higher competitive effect will have an advantage, but will still be coexistant

simDIFFcoexist=ode(y=y,times=t,func=diffSolver,parms=pDIFFcoexist)
OutputDIFFcoexist<-data.frame(time=simDIFFcoexist[,1],N1=simDIFFcoexist[,2],N2=simDIFFcoexist[,3])

ggplot(OutputDIFFcoexist,aes(x=time))+
  geom_point(aes(y=N1),color="red")+
  geom_point(aes(y=N2),color="green")+
  ylim(0,1)

Combined<-cbind.data.frame(Output1champ,Output2champ[,2:3],OutputDIFFcoexist[,2:3],OutputEQcoexist[,2:3])
names(Combined)<-c("time","Case1N1","Case1N2","Case2N1","Case2N2","Case3N1","Case3N2","Case4N1","Case4N2")
library(reshape2)
Combined<-melt(Combined,id="time")
Combined$Case<-gsub("N[0-9]","",Combined$variable)
Combined$N<-gsub("Case[0-9]","",Combined$variable)
ggplot(Combined, aes(x=time,y=value,color=Case,shape=N))+
  geom_point()


