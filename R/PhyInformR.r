##libraries
library(ape)
library(splines)
library(gplots)
library(RColorBrewer)
require(phytools)
library(foreach)
library(iterators)
library(geiger)
library(doParallel)
library(gridExtra)
library(ggplot2)
library(hexbin)
library(PBSmodelling)
## modern deps
library(dplyr)
library(purrr)
library(tibble)
library(plotly)
library(htmlwidgets)
#set the number of cores
registerDoParallel(cores=8)

##internal functions, just read these in###########################
###########################
###########################
###########################
###Need this to get individual point at a time of profile
#this function is based on Townsend 2007
site.summer<-function(rate.vector,time)
{
  length(rate.vector)->calculation.length
  at.site<-matrix(ncol=calculation.length)
  for(i in 1:calculation.length)
    
    
  {
    rate.vector[i]->current
    16*current*current*time*exp(-4*current*time)->at.site[i]
  }
  sum(at.site)->inform.at.time
  return(inform.at.time)
}

##another internal function
###########################
###########################
###########################

get.ind.sites<-function(rate.output,breaks)
{
  
  rate.output->rates
  length(rates)->vector.length
  c(1:vector.length)->numbers
  cbind(numbers,rates)->unsorted.matrix
  length(breaks[,1])->n
  length(rates)->limit
  matrix(ncol=n, nrow=limit)->extracted.sites
  matrix(ncol=n)->names.of.columns
  for(i in 1:n)
  {
    ###this part looks through the breaks and extracts the site numbers for each user specified bin
    breaks[i,]->upper.lower
    upper.lower[1]->lower
    upper.lower[2]->upper
    which(rates>=lower)->lista
    which(rates<=upper)->listb
    ####get the list of sites, which are bigger than lower bound but smaller than upper bound
    lista[(lista%in%listb)]->numbers
    length(numbers)->data.length
    limit-data.length->filler	
    rep("Na",filler)->fill
    c(numbers,fill)->output
    output-> extracted.sites[,i]
    
  }
  ###assign column names
  for(i in 1:n)
  {
    string1="Charset_"
    string2=paste(string1,i,sep="")
    string3=paste(string2,":",sep="")
    names.of.columns[,i]<-string3
  }
  colnames(extracted.sites)<-names.of.columns
  as.data.frame(extracted.sites)->ES
  return(ES)
  
}	

###internal function #3
###########################
###########################
###########################
###########################

inform.profile.generator2<-function(use.rates,tree)
{
  branching.times(tree)->btimes
  c(0,btimes)->btimes2
  sort(btimes2)->sorted.btimes
  length(btimes2)->branching.points
  length(use.rates)->calculation.length
  
  inform.at.time<-matrix(ncol=branching.points)
  for(i in 1:branching.points)
  {
    sorted.btimes[i]->btime
    site.summer(use.rates,btime)->inform.at.time[i]
    
    
  }
  inform.at.time->close
  return(close)
}
#' @export
defined.multi.profile<-function(rate.vector,tree,breaks)
{
  
  length(rate.vector)->n
  branching.times(tree)->btimes
  c(0,btimes)->btimes2
  sort(btimes2)->sorted.btimes
  length(btimes2)->branching.points
  
  
  length(breaks[,1])->n.parts
  close<-matrix(ncol=branching.points,nrow=n.parts)
  for (i in 1:n.parts)
  {
    
    min(breaks[i,]):max(breaks[i,])->part
    as.matrix(part)->partx
    partx[partx%in%1:n]->part.check
    as.numeric(part.check)->part.check
    rate.vector->rates
    rates[part.check]->part.current
    inform.profile.generator2(part.current,tree)->close[i,]
    
  }
  
  
  
  rbind(sorted.btimes,close)->closer
  return(closer)
}

#' @export
Approximator<-function(t,t0,rateVector,s)	
{	
  rateVector->rv
  currentProbability<-matrix(nrow=length(rv), ncol=1)
  Expectationxinnersum1<-c(0)
  Expectationxinnersum2<-c(0)
  Expectationy<-c(0)
  Expectationy2<-c(0)
  ExpectationX1Y<-c(0)
  ExpectationSQROOTX1Y<-c(0)
  length(rv)->n
  
  ###Loop calculations and variance
  for(i in 1:n)
  {
    rv[i,]->rateVector2
    npnl<-pnl(rateVector2,t,t0,s)
    npro<-prother(rateVector2,t,t0,s)
    npsnr<-psnr(rateVector2,t,t0,s)
    
    Expectationy<-Expectationy+npsnr
    Expectationxinnersum1<-Expectationxinnersum1+npnl
    Expectationxinnersum2<-Expectationxinnersum2+npnl*npnl
    Expectationy2<-Expectationy2+ npsnr* npsnr
    ExpectationX1Y<-ExpectationX1Y+ npsnr* npnl
    ExpectationSQROOTX1Y<-ExpectationSQROOTX1Y+ npsnr*sqrt(npnl)
  }
  
  
  Expectationx<-Expectationxinnersum1+sqrt((Expectationxinnersum1/pi))
  
  
  Expectation<- Expectationy-Expectationx
  variancey<- Expectationy-Expectationy2
  variancex<-((pi-1)/pi)*Expectationxinnersum1-Expectationxinnersum2
  variance<-variancey+variancex-2*ExpectationX1Y-(2/sqrt(pi)) * ExpectationSQROOTX1Y	
  
  rnorm(n, mean=Expectation, sd=sqrt(variance))->ndistr
  princtree<-pnorm(-0.5,mean=Expectation, sd=sqrt(variance))
  prpolytomy<-pnorm(0.5,mean=Expectation, sd=sqrt(variance))-pnorm(-0.5,mean=Expectation, sd=sqrt(variance))
  prcortree=1-pnorm(0.5,mean=Expectation, sd=sqrt(variance))
  c("Probabilty Correct", "Probability Polytomy", "Probability Incorrect" )->labels
  c(prcortree,prpolytomy,princtree)->values
  labels->names(values)
  return(values)
}

psnr<-function(lambda,t,t0,s)
{
  (-1/s^3 + 1/s^2 + (3/s^3 -1/s^2 -2/s +1 + (-4/s^2 + 4/s -1)*exp(-t0*lambda)) *exp(-(4*s)/(s-1)*t*lambda) + (-8/s^3+4/s^2  +(8/s^2-4/s) *exp (-t0*lambda)) *exp((-3*s)/(s-1)*t*lambda) +(6/s^3 -4/s^2 +2/s -4/s^2 *exp(-t0*lambda) )*exp((-2*s)/(s-1)*t*lambda) )->psnr.value
  return(psnr.value)}


pnl<-function(lambda,t,t0,s)
{pnl.value<-( -1/s^3 +1/s^2 + (3/s^3-1/s^2 + (-4/s^2+2/s) * exp(-t0*lambda) ) *exp(((-4*s)/(s-1))*t*lambda) + (-8/s^3+4/s^2+ (8/s^2-4/s) * exp(-t0*lambda)) * exp((-3*s)/(s-1)*t*lambda) + (6/s^3 - 4/s^2 + (-4/s^2+2/s) * exp(-t0*lambda))*exp((-2*s/(s-1))*t*lambda) )
return(pnl.value)
}

pnL2<-function(lambda,t,t0,s)
{pnl.value<-(-1/s^3+1/s^2+ (3/s^3-1/s^2+ (-4/s^2+2/s) *exp(-t0*lambda) ) *exp(((-4*s)/(s-1))*t*lambda) + (-8/s^3+4/s^2+(8/s^2-4/s)*exp(-t0*lambda)) * exp((-3*s)/(s-1)*t*lambda) +(6/s^3 - 4/s^2 + (-4/s^2+2/s) *exp(-t0*lambda))*exp((-2*s/(s-1))*t*lambda) )
return(pnl.value)
}

prother<-function(lambda,t,t0,s)
{prother.value<-1-pnL2(lambda,t,t0,s)-pnl(lambda,t,t0,s)-psnr(lambda,t,t0,s)
return(prother.value)}

pnl2<-function(lambda,t,t0,s)
{pnl.value<-(-1/s^3+1/s^2+ (3/s^3-1/s^2+ (-4/s^2+2/s) *exp(-t0*lambda) ) *exp(((-4*s)/(s-1))*t*lambda) + (-8/s^3+4/s^2+(8/s^2-4/s)*exp(-t0*lambda)) * exp((-3*s)/(s-1)*t*lambda) +(6/s^3 - 4/s^2 + (-4/s^2+2/s) *exp(-t0*lambda))*exp((-2*s/(s-1))*t*lambda) )
return(pnl.value)
}

prother<-function(lambda,t,t0,s)
{prother.value<-1-pnL2(lambda,t,t0,s)-pnl(lambda,t,t0,s)-psnr(lambda,t,t0,s)
return(prother.value)}
Approximator.lite<-function(t,t0,rateVector,s)	
{	
  rateVector->rv
  currentProbability<-matrix(nrow=length(rv), ncol=1)
  Expectationxinnersum1<-c(0)
  Expectationxinnersum2<-c(0)
  Expectationy<-c(0)
  Expectationy2<-c(0)
  ExpectationX1Y<-c(0)
  ExpectationSQROOTX1Y<-c(0)
  length(rv)->n
  
  ###Loop calculations and variance
  for(i in 1:n)
  {
    rv[i,]->rateVector2
    npnl<-pnl(rateVector2,t,t0,s)
    npro<-prother(rateVector2,t,t0,s)
    npsnr<-psnr(rateVector2,t,t0,s)
    
    Expectationy<-Expectationy+npsnr
    Expectationxinnersum1<-Expectationxinnersum1+npnl
    Expectationxinnersum2<-Expectationxinnersum2+npnl*npnl
    Expectationy2<-Expectationy2+ npsnr* npsnr
    ExpectationX1Y<-ExpectationX1Y+ npsnr* npnl
    ExpectationSQROOTX1Y<-ExpectationSQROOTX1Y+ npsnr*sqrt(npnl)
  }
  
  
  Expectationx<-Expectationxinnersum1+sqrt((Expectationxinnersum1/pi))
  
  
  Expectation<- Expectationy-Expectationx
  variancey<- Expectationy-Expectationy2
  variancex<-((pi-1)/pi)*Expectationxinnersum1-Expectationxinnersum2
  variance<-variancey+variancex-2*ExpectationX1Y-(2/sqrt(pi)) * ExpectationSQROOTX1Y	
  
  rnorm(n, mean=Expectation, sd=sqrt(variance))->ndistr
  princtree<-pnorm(-0.5,mean=Expectation, sd=sqrt(variance))
  prpolytomy<-pnorm(0.5,mean=Expectation, sd=sqrt(variance))-pnorm(-0.5,mean=Expectation, sd=sqrt(variance))
  prcortree=1-pnorm(0.5,mean=Expectation, sd=sqrt(variance))
  c("Probabilty Correct", "Probability Polytomy", "Probability Incorrect" )->labels
  c(prcortree,prpolytomy,princtree)->values
  labels->names(values)
  return(prcortree)
}

#' @export
space.maker<-function(rateVector,t,s)
{
  t/20->by.this
  seq(by.this,t-0.0001,by=by.this)->lilts
  rowspace<-matrix(nrow=1,ncol=length(lilts))
  for (i in 1:length(lilts))
  {
    lilts[i]->to
    Approximator.lite(t,to,rateVector,s)->rowspace[i]
  }
  return(rowspace)
}

#' @export
space.maker.narrow<-function(rateVector,t,s)
{
  t/2->halft
  halft/20->by.this
  seq(by.this, halft-0.0001,by=by.this)->lilts
  rowspace<-matrix(nrow=1,ncol=length(lilts))
  for (i in 1:length(lilts))
  {
    lilts[i]->to
    Approximator.lite(t,to,rateVector,s)->rowspace[i]
  }
  return(rowspace)
}




##generates informativeness output like phydesign
inform.profile.generator<-function(rate.vector,tree)
{
  
  branching.times(tree)->btimes
  c(0,btimes)->btimes2
  sort(btimes2)->sorted.btimes
  length(btimes2)->branching.points
  length(rate.vector)->calculation.length
  
  inform.at.time<-matrix(ncol=branching.points)
  for(i in 1:branching.points)
  {
    sorted.btimes[i]->btime
    site.summer(rate.vector,btime)->inform.at.time[i]
    
    
  }
  rbind(sorted.btimes,inform.at.time)->close
  return(close)
}

####This part will get all the points with the rate vector already computed from other functions
#' @export
informativeness.profile<-function(rate.vector, tree, codon="FALSE", values="display")
{
  
  branching.times(tree)->btimes
  c(0,btimes)->btimes2
  
  
  if (codon=="FALSE"){
    inform.profile.generator(rate.vector,tree)->close
    
    
    close[1,]->sorted.btimes
    close[2,]->inform.at.time
    round(max(btimes))->upper
    upper/5->by.this
    round(max(inform.at.time),digits=2)->uppery
    uppery/10->by.y
    
    yy <-predict(interpSpline(sorted.btimes, inform.at.time))
    mat<- matrix(c(1:2),nrow=2,ncol=1)
    graphics::layout(mat=mat,heights=c(250,300))
    graphics::par(mar=c(0,0,0,0), oma=c(5,5,1,1))
    #par(bg = "white")   
    #split.screen(c(2,1))
    #screen(1)
    #plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
    ##coord are left,right,bottom,top from 0 to 1
    
    graphics::par(plt=c(0,0.9,0.2,0.99))
    plot(tree,show.tip.label=FALSE,direction="l")
    ####Lower corner, note that the pi is offset to mirror the trees end
    
    graphics::par(plt=c(0.027,0.9,0,0.99))
    
    plot(sorted.btimes,inform.at.time,pch=NA_integer_,axes=FALSE, ylim=c(0,uppery+(uppery*.15)), xlim=c(0,upper))
    axis(1, at = seq(0, upper, by = by.this), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
    axis(2, at = seq(0, uppery, by = by.y), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
    
    
    lines(yy, pch=NA_integer_, col="blue",lty=1,)
    legend("topright",y=NULL,c("PI of Locus"),lty=1,col="blue",lwd=2,title="PI Profile")
    #plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
    #close.screen(all = TRUE)
    
  }
  if (codon=="TRUE")
  {
    pos1 <- rate.vector[seq(1, length(rate.vector), 3)]
    pos2 <- rate.vector[seq(2, length(rate.vector), 3)]
    pos3 <- rate.vector[seq(3, length(rate.vector), 3)]
    length(btimes2)->branching.points
    
    close2<-matrix(ncol=branching.points,nrow=3)
    
    inform.profile.generator2(pos1,tree)->close2[1,]
    inform.profile.generator2(pos2,tree)->close2[2,]
    inform.profile.generator2(pos3,tree)->close2[3,]
    
    
    inform.profile.generator(rate.vector,tree)->close
    
    
    close[1,]->sorted.btimes
    close[2,]->inform.at.time
    round(max(btimes))->upper
    upper/5->by.this
    sort(btimes2)->sortedbtimes2
    round(max(close2),digits=2)->uppery
    uppery/10->by.y
    
    mat<- matrix(c(1:2),nrow=2,ncol=1)
    graphics::layout(mat=mat,heights=c(250,300))
    graphics::par(mar=c(0,0,0,0), oma=c(5,5,1,1))
    #par(bg = "white")   
    #split.screen(c(2,1))
    #screen(1)
    #plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
    ##coord are left,right,bottom,top from 0 to 1
    
    graphics::par(plt=c(0,0.9,0.2,0.99))
    plot(tree,show.tip.label=FALSE,direction="l")
    ####Lower corner, note that the pi is offset to mirror the trees end
    
    graphics::par(plt=c(0.027,0.9,0,0.99))
    
    plot(sorted.btimes,close2[3,],pch=NA_integer_,axes=FALSE, ylim=c(0,uppery), xlim=c(0,upper))
    
    axis(1, at = seq(0, upper, by = by.this), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),ylab="Time from Present")
    axis(2, at = seq(0, uppery, by = by.y), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),  xlab="Phylogenetic Informativeness")
    c("black","blue","gray","green","purple","brown","azure","red","yellow")->colors
    c("part1","part2","part3","part4","part5","part6","part7","part8","part9")->leglab
    c(1,2,3,1,2,3,1,2,3)->style
    c(2,2,3,2,2,3,2,2,3)->thickness
    legend("topright",y=NULL,leglab[1:length(close2[,1])],lty=style,col=colors,lwd=thickness,title="Partition PI Profile")
    for (i in 1:length(close2[,1])){
      close2[i,]->inform.at.time	
      yy <-predict(interpSpline(sortedbtimes2, inform.at.time))
      
      lines(yy, pch=NA_integer_, col=colors[i],lty=style[i],)
      
    }
    resetGraph(reset.mf=TRUE)
    rbind(sorted.btimes,close2)->closer
    if (values=="display"){
      return(closer)} else if (values=="off"){
        return("done")
      }
  }
  #return(close)
}

####For user defined informativeness profiles, note that this has a maximum limit of X since the plot will become unreadable
#' @export
multi.profile<-function(rate.vector,tree,breaks,values="display")
{
  length(rate.vector)->n
  branching.times(tree)->btimes
  c(0,btimes)->btimes2
  sort(btimes2)->sorted.btimes
  length(btimes2)->branching.points
  
  get.ind.sites(rate.vector,breaks)->ES
  length(breaks[,1])->n.parts
  close<-matrix(ncol=branching.points,nrow=n.parts)
  for (i in 1:n.parts)
  {
    
    ES[,i]->part
    as.matrix(part)->partx
    partx[partx%in%1:n]->part.check
    as.numeric(part.check)->part.check
    rate.vector->rates
    rates[part.check]->part.current
    inform.profile.generator2(part.current,tree)->close[i,]
    
  }
  
  
  ####now draw the profile######################
  ########First set the x and y axis bounds##########
  round(max(btimes))->upper
  upper/5->by.this
  
  round(max(close))->uppery
  uppery/10->by.y
  
  
  
  mat<- matrix(c(1:2),nrow=2,ncol=1)
  graphics::layout(mat=mat,heights=c(250,300))
  graphics::par(mar=c(0,0,0,0), oma=c(5,5,1,1))
  #par(bg = "white")   
  #split.screen(c(2,1))
  #screen(1)
  #plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
  ##coord are left,right,bottom,top from 0 to 1
  
  graphics::par(plt=c(0,0.9,0.2,0.99))
  plot(tree,show.tip.label=FALSE,direction="l")
  ####Lower corner, note that the pi is offset to mirror the trees end
  
  graphics::par(plt=c(0.027,0.9,0,0.99))
  
  plot(sorted.btimes,close[1,],pch=NA_integer_,axes=FALSE, ylim=c(0,uppery), xlim=c(0,upper))
  
  axis(1, at = seq(0, upper, by = by.this), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),ylab="Time from Present")
  axis(2, at = seq(0, uppery, by = by.y), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),  xlab="Phylogenetic Informativeness")
  c("black","blue","gray","green","purple","brown","azure","red","yellow")->colors
  c("part1","part2","part3","part4","part5","part6","part7","part8","part9")->leglab
  c(1,2,3,1,2,3,1,2,3)->style
  c(2,2,3,2,2,3,2,2,3)->thickness
  legend("topright",y=NULL,leglab[1:n.parts],lty=style,col=colors,lwd=thickness,title="Partition PI Profile")
  for (i in 1:n.parts){
    close[i,]->inform.at.time	
    yy <-predict(interpSpline(sorted.btimes, inform.at.time))
    
    lines(yy, pch=NA_integer_, col=colors[i],lty=style[i],)
    
  }
  #resetGraph(reset.mf=TRUE)
  rbind(sorted.btimes,close)->closer
  if (values=="display"){
    return(closer)} else {
      print("done")
    }
}


###this gets the output of all three positions
#c("times","pos1","pos2","pos3")->rownames

#rbind(close,inform.at.time2)->closer
#rbind(closer,inform.at.time3)->cLoser
#row.names(cLoser)<-rownames
#return(cLoser)



#}	


####For user defined informativeness profiles, note that this has a maximum limit of X since the plot will become unreadable
#' @export
defined.multi.profile<-function(rate.vector,tree,breaks,values="display")
{
  length(rate.vector)->n
  branching.times(tree)->btimes
  c(0,btimes)->btimes2
  sort(btimes2)->sorted.btimes
  length(btimes2)->branching.points
  
  
  length(breaks[,1])->n.parts
  close<-matrix(ncol=branching.points,nrow=n.parts)
  for (i in 1:n.parts)
  {
    
    min(breaks[i,]):max(breaks[i,])->part
    as.matrix(part)->partx
    partx[partx%in%1:n]->part.check
    as.numeric(part.check)->part.check
    rate.vector->rates
    rates[part.check]->part.current
    inform.profile.generator2(part.current,tree)->close[i,]
    
  }
  
  
  ####now draw the profile######################
  ########First set the x and y axis bounds##########
  round(max(btimes))->upper
  upper/5->by.this
  
  round(max(close))->uppery
  uppery/10->by.y
  
  
  
  mat<- matrix(c(1:2),nrow=2,ncol=1)
  graphics::layout(mat=mat,heights=c(250,300))
  graphics::par(mar=c(0,0,0,0), oma=c(5,5,1,1))
  #par(bg = "white")   
  #split.screen(c(2,1))
  #screen(1)
  #plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
  ##coord are left,right,bottom,top from 0 to 1
  
  graphics::par(plt=c(0,0.9,0.2,0.99))
  plot(tree,show.tip.label=FALSE,direction="l")
  ####Lower corner, note that the pi is offset to mirror the trees end
  
  graphics::par(plt=c(0.027,0.9,0,0.99))
  
  plot(sorted.btimes,close[1,],pch=NA_integer_,axes=FALSE, ylim=c(0,uppery), xlim=c(0,upper))
  
  axis(1, at = seq(0, upper, by = by.this), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),ylab="Time from Present")
  axis(2, at = seq(0, uppery, by = by.y), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0),  xlab="Phylogenetic Informativeness")
  c("black","blue","gray","green","purple","brown","azure","red","yellow")->colors
  c("part1","part2","part3","part4","part5","part6","part7","part8","part9")->leglab
  c(1,2,3,1,2,3,1,2,3)->style
  c(2,2,3,2,2,3,2,2,3)->thickness
  legend("topright",y=NULL,leglab[1:n.parts],lty=style,col=colors,lwd=thickness,title="Partition PI Profile")
  for (i in 1:n.parts){
    close[i,]->inform.at.time	
    yy <-predict(interpSpline(sorted.btimes, inform.at.time))
    
    lines(yy, pch=NA_integer_, col=colors[i],lty=style[i],)
    
  }
  #resetGraph(reset.mf=TRUE)
  rbind(sorted.btimes,close)->closer
  if (values=="display"){
    return(closer)} else {
      print("done")
      
    }
}
NodeWalker<-function(tree)
{
  
  #first figure out which tree$edge rows are just internal nodes
  rows <- which(tree$edge[,2]>length(tree$tip.label))
  #tree$edge.length[rows]
  #this gets the internal node labels (parent and daughter) and internode length
  dist <- cbind(tree$edge[rows,],tree$edge.length[rows])
  #gets the branchingtimes of parent node by matching column 1
  parent_times <- branching.times(tree)[match(dist[,1],names(branching.times(tree)))]
  #gets the branchingtimes of daughter node by matching column 1
  daughter_times	 <- branching.times(tree)[match(dist[,2],names(branching.times(tree)))]
  dist <- cbind(dist,parent_times,daughter_times)
  colnames(dist) <- c("parent_node","daughter_node","edge_length","p_node_time","d_node_time")		
  return(dist)
}

###Plot the space
#' @export
PlotTreeSI<-function(tree,ratevector,s)
{{
  #First get x axis
  NodeWalker(tree)->nodes
  nodes[,4]->parentn
  nodes[,5]->daughtern
  probs<-matrix(ncol=3,nrow=length(parentn))
  for (i in 1:length(parentn))
  {
    parentn[i]->parentValue
    daughtern[i]->t
    parentValue-t->t0
    Approximator(t,t0,ratevector,s)->probs[i,]
    
  }		
  probs[,1]->correct
  probs[,3]->incorrect
  
  mat<- matrix(c(1:2),nrow=2,ncol=1)
  graphics::layout(mat=mat,heights=c(250,300))
  graphics::par(mar=c(0,0,0,0), oma=c(5,5,1,1))
  graphics::par(plt=c(0,0.9,0.2,0.99))
  plot(tree,show.tip.label=FALSE,direction="l")
  graphics::par(plt=c(0.027,0.9,0,0.99))
  
  ###upperx is the xlim in case you want to zoom in or expand for another rate vector
  plot(parentn,correct, xlim=c(0, max(parentn)),ylim=c(0,1), col="white", pch=17)
  #points(daughtern,correct, col="blue", pch=17)
  segments(parentn,correct,daughtern,correct, col="blue")
  
  #points(parentn,incorrect, bg=312, pch=25)
  #points(daughtern,incorrect, bg=312, pch=25)
  #segments(parentn,incorrect,daughtern,incorrect, col=312)
  return(probs)
}
  resetGraph(reset.mf=TRUE)
}
#' @export
Plot.Another.TreeSI<-function(tree,ratevector,s,color,type)
{
  #First get x axis
  NodeWalker(tree)->nodes
  nodes[,4]->parentn
  nodes[,5]->daughtern
  probs<-matrix(ncol=3,nrow=length(parentn))
  for (i in 1:length(parentn))
  {
    parentn[i]->parentValue
    daughtern[i]->t
    parentValue-t->t0
    Approximator(t,t0,ratevector,s)->probs[i,]
    
  }		
  probs[,1]->correct
  probs[,3]->incorrect
  #points(daughtern,correct, col=color, pch=17)
  #points(parentn, correct, bg=312, pch=25)
  graphics::par(plt=c(0.027,0.9,0,0.99))
  
  segments(parentn,correct,daughtern,correct, col=color,lty=type)
  
  #points(parentn,incorrect, bg=312, pch=25)
  #points(daughtern,incorrect, bg=312, pch=25)
  #segments(parentn,incorrect,daughtern,incorrect, col=312)
  return(probs)
  resetGraph(reset.mf=TRUE)
  
}

PlotManyTreeSN<-function(contree,trees,ratevector,s)
{
  
  sapply(trees,NodeWalker)->nodes
  
  
}

####uses rate output from rate.by.site function and a matrix of breaks to return the site numbers for the break matrix should be in the format lower bound in column one, upper bound in column 2 this can have as many partitions as the user wants. Best idea would be to use the hist(rate.vector[1,]) function to see the frequency distribution of site patterns to design the partitioning strategy	

####Note that this function returns site patterns as a data.frame. Need another function to deal with this output for actual nexus output

get.ind.sites<-function(rate.output,breaks)
{
  
  rate.output->rates
  length(rates)->vector.length
  c(1:vector.length)->numbers
  cbind(numbers,rates)->unsorted.matrix
  length(breaks[,1])->n
  length(rates)->limit
  matrix(ncol=n, nrow=limit)->extracted.sites
  matrix(ncol=n)->names.of.columns
  for(i in 1:n)
  {
    ###this part looks through the breaks and extracts the site numbers for each user specified bin
    breaks[i,]->upper.lower
    upper.lower[1]->lower
    upper.lower[2]->upper
    which(rates>=lower)->lista
    which(rates<=upper)->listb
    ####get the list of sites, which are bigger than lower bound but smaller than upper bound
    lista[(lista%in%listb)]->numbers
    length(numbers)->data.length
    limit-data.length->filler	
    rep("Na",filler)->fill
    c(numbers,fill)->output
    output-> extracted.sites[,i]
    
  }
  ###assign column names
  for(i in 1:n)
  {
    string1="Charset_"
    string2=paste(string1,i,sep="")
    string3=paste(string2,":",sep="")
    names.of.columns[,i]<-string3
  }
  colnames(extracted.sites)<-names.of.columns
  as.data.frame(extracted.sites)->ES
  return(ES)
  
}	

#to run and log output to cluster, note image name and filename cannot be identical!!!!
###Sample Input: cluster.signal.noise(86, 91, cytBM, 10, filename="tarsius", imagename="tarsius.pdf")
#' @export	
cluster.signal.noise<-function(t, t0, rateVector, nsims,s, filename,imagename, image="FALSE")
{
  signal.noise.multimix(t,t0,rateVector, nsims,s)->currentprobdist
  normdprobdist<-(currentprobdist/nsims)
  ###probability of yielding correct tree
  length(rateVector)->n
  2*n+1->max.bound
  n+2->start.right
  n+1->poly
  ###probability of getting the right tree
  normdprobdist[start.right:max.bound]->right.signal
  sum(right.signal)->right.signal
  #return(right.signal)
  ###probability of polytomy
  normdprobdist[poly]->polytomy
  #return(polytomy)
  ###probability of wrong tree
  normdprobdist[1:n]->wrong.tree
  sum(wrong.tree)->false.knowledge
  #return(false.knowledge)
  ###odds ratio of correct vs. incorrect tree
  right.signal/false.knowledge->odds.of.recovery
  
  ###odds of correct vs incorrect OR polytomy
  wrong.tree+polytomy->bogus
  right.signal/bogus->odds.of.anything
  
  ###odds of no polytomy
  1-polytomy->odds.resolving
  ###plot using function below
  graph.signal.noise(currentprobdist, rateVector, imagename, image)
  
  ##return values
  sig.noise<-cbind(right.signal,polytomy,false.knowledge,odds.of.recovery,odds.of.anything,odds.resolving)
  colnames(sig.noise)<-c("P_correct", "P_poly","P_wrong", "odds_correctvswrong", "odds_corrvswrong/poly", "odds_resolving")
  write.table(sig.noise[1,],  file=filename)
}


parallel.multimixA<-function (t,t0, ratevector, nsims,s){
  ratevector->rv
  n<-length(rv)
  matrix(nrow=n*2+1, ncol=nsims)->shrill.mess
  foreach(q=1:nsims, .combine=cbind)%dopar%
    {
      CurrentProbabilityDistribution(ratevector, t, t0,s)->shrill.mess[,q]
      
      
    }	
  
}

parallel.multimixfull<-function(t,t0,ratevector, nsims, s)
{
  parallel.multimixA(t,t0,ratevector, nsims, s)->temp
  rowSums(temp)->currentprobdistro
  return(currentprobdistro)
  
  
}
#to run and log output to cluster, note image name and filename cannot be identical!!!!
#' @export
parallel.cluster.signal.noise<-function(t, t0, rateVector, nsims,s, filename,imagename, image="TRUE")
{
  parallel.multimixfull(t,t0,rateVector, nsims,s)->currentprobdist
  normdprobdist<-(currentprobdist/nsims)
  ###probability of yielding correct tree
  length(rateVector)->n
  2*n+1->max.bound
  n+2->start.right
  n+1->poly
  ###probability of getting the right tree
  normdprobdist[start.right:max.bound]->right.signal
  sum(right.signal)->right.signal
  #return(right.signal)
  ###probability of polytomy
  normdprobdist[poly]->polytomy
  #return(polytomy)
  ###probability of wrong tree
  normdprobdist[1:n]->wrong.tree
  sum(wrong.tree)->false.knowledge
  #return(false.knowledge)
  ###odds ratio of correct vs. incorrect tree
  right.signal/false.knowledge->odds.of.recovery
  
  ###odds of correct vs incorrect OR polytomy
  wrong.tree+polytomy->bogus
  right.signal/bogus->odds.of.anything
  
  ###odds of no polytomy
  1-polytomy->odds.resolving
  ###plot using function below
  graph.signal.noise(currentprobdist, rateVector, imagename, image)
  
  ##return values
  sig.noise<-cbind(right.signal,polytomy,false.knowledge,odds.of.recovery,odds.of.anything,odds.resolving)
  colnames(sig.noise)<-c("P_correct", "P_poly","P_wrong", "odds_correctvswrong", "odds_corrvswrong/poly", "odds_resolving")
  write.table(sig.noise[1,],  file=filename)
  
  if (image=="TRUE"){
    return(sig.noise[1,])
  } else if (image=="FALSE"){
    return("done")
  }
}

###Solve for x you need for current probabilty distribution	
##s is the number of states	
ExMaker<-function(t,t0,rateVector,s)
{	
  rateVector->rv
  currentProbability<-matrix(nrow=length(rv), ncol=1)
  nwone<-matrix(nrow=length(rv), ncol=1)
  nwtwo<-matrix(nrow=length(rv), ncol=1)
  snr<-matrix(nrow=length(rv), ncol=1)
  length(rv)->n
  for(i in 1:n)
  {rv[i,]->lambda
    ###evaluate lambda	
    npro<-prother(lambda,t,t0,s)
    npsnr<-psnr(lambda,t,t0,s)
    npnl<-pnl(lambda,t,t0,s)
    npnl2<-pnl2(lambda,t,t0,s)
    ###empirical test Block
    #npro<-prother(.21,.9,.02)
    #	npsnr<-psnr(.21,.9,.02)
    #	npnl<-pnl(.21,.9,.02)
    #	npnl2<-pnl2(.21,.9,.02)
    
    ##draw random Number
    runif(1, min=0, max=1)->randomNumber
    ###create columns of values
    ###first=null
    ###testloop
    #Returner<-function(randomNumber, npsnr,npnl,npnl2,npro,matrix){	
    #n<-length(matrix)
    #	currentProbability<-matrix(nrow=length(matrix), ncol=1)
    #nwone<-0
    #nwtwo<-0
    #snr<-0
    if (randomNumber<npro){
      next}else
        ###second=signal noise ratio
        if(randomNumber<npro+npsnr){1->snr[i,]
        }else
          ###third=nwone
          if(randomNumber<npro+npsnr+npnl){1->nwone[i,]
          }else			
            ###fourth=nwtwo
            if(randomNumber<npro+npsnr+npnl){1->nwtwo[i,] 
            }
  }
  na.omit(snr)->snr
  sum(snr)->snr
  na.omit(nwone)->nwone
  sum(nwone)->nwone
  na.omit(nwtwo)->nwtwo
  sum(nwtwo)->nwtwo
  cbind(nwone, nwtwo)->wrong
  max(wrong)->wronger	
  currentProbability2<-n+1+snr-wronger
  return(currentProbability2)}

CurrentProbabilityDistribution<-function(rateVector, t, t0,s)
{
  rateVector->rv
  n<-length(rv)
  ExMaker(t,t0, rateVector,s)->x
  ####create zero vector
  mat.or.vec(n*2+1,1)->zero.table
  as.matrix(zero.table)->zero.table.mat
  zero.table.mat[x,1]<-1
  return(zero.table.mat)
}

###Put all into one matrix
signal.noise.multimix<-function(t,t0,rateVector, nsims,s)
{
  rateVector->rv
  n<-length(rv)
  matrix(nrow=n*2+1, ncol=nsims)->shrill.mess
  for(i in 1:nsims)
  {
    CurrentProbabilityDistribution(rateVector, t, t0,s)->shrill.mess[,i]
  }
  rowSums(shrill.mess)->current.prob.dist
  return(current.prob.dist)}

#internal drawing function for Terminal Runs

graph.signal.noise<-function(currentprobdist, rateVector, filename, image="TRUE")
{
  currentprobdist->matrix.of.noise
  length(rateVector)->n
  2*n+1->max.bound
  n+2->start.right
  n+1->poly
  ###isolate_each
  
  matrix.of.noise[poly]->polytomy
  
  ###isolate_Non_zero_grey
  matrix.of.noise[start.right:max.bound]->green.side
  green.side!=0->nozero1
  green.side[nozero1]->blue.plot
  ###isolate_Non_zero_blue
  matrix.of.noise[1:n]->wrong.side
  wrong.side!=0->nozero2
  wrong.side[nozero2]->grey.plot
  
  ###the x y for all
  c(grey.plot,polytomy,blue.plot)->plotter
  1:length(plotter)->plotterx
  
  ##the wrong side and poly parts
  length(grey.plot)->wrong.side.length
  1:wrong.side.length->wrong.side.x
  1+wrong.side.length->poly.location
  
  ###the prob of correct signal parts
  as.numeric(max(plotterx))->upperbound
  1+poly.location->bluestart
  bluestart:upperbound->correct.x
  
  
  if (image=="TRUE"){
    plot(plotterx,plotter, type="h", lwd=6, xaxt="n", bty="l", ylab="Frequency", xlab="Signal Noise Plot")
    if(length(polytomy)>=1)
    {
      lines(poly.location, polytomy, type="h", col="black", xaxt="n", lwd=12)}
    if(length(grey.plot)>=1){
      lines(wrong.side.x, grey.plot, type="h", col="grey", xaxt="n", lwd=12)}
    if(length(blue.plot)>=1)
    {
      lines(correct.x, blue.plot, type="h", col="blue", xaxt="n", lwd=12)}
  } else if (image=="FALSE") {
    
    
    
    ##disregard
    #bin<-as.numeric(cut(polytomy,combo.time$breaks))
    #plot(combo.time, col=replace(rep("white", length(combo.time$breaks)-1), bin, "blue"))
    
    pdf(file=filename, height=8,width=8)
    
    plot(plotterx,plotter, type="h", lwd=6, xaxt="n", bty="l", ylab="Frequency", xlab="Signal Noise Plot")
    if(length(polytomy)>=1)
    {
      lines(poly.location, polytomy, type="h", col="black", xaxt="n", lwd=12)}
    if(length(grey.plot)>=1){
      lines(wrong.side.x, grey.plot, type="h", col="grey", xaxt="n", lwd=12)}
    if(length(blue.plot)>=1)
    {
      lines(correct.x, blue.plot, type="h", col="blue", xaxt="n", lwd=12)}
    dev.off()
  }}
#' @export
allmodel.signal.noise<-function(a,b,c,d,e,f,internode,Pi_T,Pi_C,Pi_A,Pi_G, rate_vector)
{
  rate_vector->rr
  
  ##Legacy call to Su_et_al.py. 
  
  #print(internode)
  #paste("-i","--internode",sep=" ")->inttemp
  #paste(inttemp,internode[1],sep=" ")->int1
  #paste("",rr, sep=" ")->rrr
  #as.numeric(rrr)->ra
  #paste("-r","--ratevector", sep=" ")->rtemp
  #paste(rtemp ,ra[1], sep= " ")->rrrr
  #paste("",ra[2:length(ra)], sep= "")->rara
  #paste("python", "./Su_et_al.py", sep=" ")->start
  #as.vector(c(a,b,c,d,e,f, Pi_T,Pi_C,Pi_A,Pi_G, int1,internode[2],internode[3],internode[4],internode[5],rrrr,rara, sep=" "))->command
  #c(start, command)->go
  #paste(go, sep=" ", collapse=" ")->go2
  #system(go2, intern=TRUE)
  
  default_rate_vector<-c(0.003108, 0, 0, 0.015862, 0.000426, 0, 0.005114, 0, 0, 0.00778, 0, 0, 0.001352, 0.000862, 0, 0.000862, 0, 0, 0.001338, 0, 0, 0, 0, 0, 0.005988, 0, 0, 0.001909, 0, 0, 0.000428, 0, 0, 0.000425, 0, 0, 0, 0, 0, 0.004616, 0, 0, 0.000869, 0.000426, 0, 0.000426, 0, 0, 0.004359, 0, 0, 0.001804, 0, 0.000424, 0.003546, 0, 0, 0.003128, 0, 0.000426, 0.015036, 0, 0, 0.005353, 0, 0, 0.004196, 0, 0, 0.002492, 0, 0, 0.005289, 0, 0, 0.004892, 0, 0, 0.005636, 0, 0, 0.002417, 0, 0, 0.003462, 0, 0, 0.001851, 0.000423, 0, 0.000423, 0, 0, 0, 0.000426, 0, 0.004764, 0.001354, 0, 0.00255, 0, 0, 0.004156, 0.001327, 0, 0.004163, 0.001361, 0, 0.001852, 0, 0, 0.000907, 0, 0, 0.002612, 0, 0, 0.001807, 0, 0, 0.003801, 0, 0, 0.001846, 0, 0, 0, 0.000425, 0, 0.004002, 0, 0, 0.009891, 0.000923, 0.001814, 0.002406, 0, 0, 0.000426, 0, 0, 0, 0, 0, 0.00608, 0, 0, 0.001862, 0, 0, 0, 0, 0, 0.003463, 0, 0, 0.003651, 0, 0, 0.000426, 0, 0, 0, 0, 0, 0.003501, 0, 0, 0.000871, 0, 0, 0.005557, 0, 0, 0.001893, 0, 0, 0.000866, 0, 0, 0.001412, 0, 0, 0.004276, 0, 0, 0.002342, 0.001837, 0, 0.001837, 0, 0, 0.001802, 0, 0, 0.005248, 0, 0, 0.001803, 0, 0, 0.001822, 0, 0, 0.001918, 0, 0, 0.003073, 0, 0, 0.001323, 0, 0, 0.009164, 0.002059, 0, 0.002637, 0.000423, 0, 0.002336, 0, 0, 0.003847, 0, 0, 0.004949, 0, 0, 0.002946, 0, 0, 0.001807, 0.000426, 0.000426, 0.000426, 0, 0, 0.003, 0, 0, 0.005688,0, 0, 0.004278, 0.001811, 0.002346, 0.012034, 0, 0, 0.001409, 0.000865, 0, 0.000865, 0, 0, 0.001374, 0, 0, 0.002942, 0, 0, 0, 0, 0.000428, 0.008127, 0, 0, 0.001892, 0, 0, 0.003498, 0.001856, 0, 0.000428, 0, 0, 0.004151, 0, 0, 0.003209, 0, 0, 0.004108, 0, 0, 0.000951, 0, 0, 0.001352, 0, 0, 0.002333, 0, 0, 0.002329, 0, 0, 0.010802, 0, 0, 0.001418, 0, 0, 0.001322, 0, 0, 0.003694, 0.001999, 0.001999, 0.005564, 0, 0, 0.007526, 0, 0, 0.003692, 0, 0, 0.003083, 0.000426, 0, 0.008106, 0.001333, 0.000425, 0.003509, 0, 0, 0.009753,0.001374, 0, 0.006182, 0, 0, 0.001363, 0.000426, 0, 0.00542, 0.001324, 0.001324, 0.004788, 0, 0, 0.000428, 0, 0, 0, 0, 0, 0.006989, 0, 0, 0.006022, 0, 0, 0, 0, 0, 0.004086, 0, 0, 0.003316, 0.000423, 0, 0.003664, 0, 0, 0.005446, 0, 0, 0.005158, 0, 0, 0, 0, 0, 0.002359, 0, 0, 0, 0, 0, 0.002336, 0, 0, 0.003833, 0, 0, 0, 0.000423, 0, 0.002407, 0, 0, 0.003585, 0, 0, 0.002905, 0, 0, 0.005398, 0, 0, 0.001824, 0.000877, 0, 0.005099, 0, 0, 0.000423, 0, 0.000425, 0, 0.000894, 0, 0.010747, 0.002316, 0, 0.005676, 0.000428,0, 0.004035, 0, 0, 0.003574, 0, 0.001347, 0.00183, 0, 0, 0.00385, 0.000876, 0.000424, 0.001835, 0.000428, 0, 0.000428, 0.001382, 0, 0.005137, 0.000423, 0, 0.003118, 0.00087, 0, 0.003728, 0, 0, 0.00405, 0, 0, 0.00087, 0.00134, 0, 0.00134, 0, 0, 0.004536, 0.000425, 0, 0.002412, 0, 0, 0.007825, 0.000874, 0.000424, 0.001347,0.001857, 0, 0.000878, 0.001349, 0, 0.002333, 0.000426, 0.000426, 0.000426, 0.000426, 0.000426, 0.002371, 0, 0, 0.00296, 0.001823, 0, 0.002912, 0.001813, 0.000884, 0.007372, 0, 0, 0.002954, 0.001373, 0,0.001893, 0, 0, 0.001343, 0, 0, 0, 0.003074, 0.000426, 0.006498, 0,0, 0, 0.001854, 0, 0.007631, 0, 0, 0.003719, 0.000874, 0.000426, 0.005504, 0, 0, 0.004131, 0, 0, 0.003597, 0.000869, 0.000428, 0.001836, 0, 0, 0.000423, 0, 0, 0.005095, 0, 0, 0.008057, 0.001426, 0, 0.001426, 0, 0, 0.003061, 0, 0, 0.00459, 0, 0, 0.004175, 0, 0, 0.005326, 0, 0.000441, 0.004984, 0.002378, 0, 0.003662, 0.000428, 0.000435, 0.000423, 0, 0, 0.002364, 0, 0, 0.004486, 0, 0, 0.003928, 0.000871, 0, 0.00354, 0, 0, 0.00568, 0, 0.000874, 0.004861, 0, 0, 0.00283, 0, 0, 0.001818, 0, 0, 0.004818, 0.001412, 0, 0.000881, 0, 0, 0.000423, 0, 0.000426, 0.006653, 0.001377, 0, 0.007102, 0.001848,0, 0.003496, 0, 0, 0.000423, 0, 0, 0.005761, 0, 0, 0.006607, 0.000431, 0, 0.009237, 0.000425, 0, 0.004134, 0, 0, 0.003539, 0, 0, 0.004863, 0, 0, 0.006153, 0, 0, 0.001959, 0, 0, 0.000884, 0.000423, 0, 0.005317, 0, 0, 0.002122, 0, 0, 0, 0, 0, 0.001811, 0.000426, 0.000426, 0.000428, 0, 0, 0.00602, 0, 0, 0.002454, 0, 0, 0.003476, 0, 0, 0.004903, 0, 0, 0.000428, 0, 0, 0.001404, 0, 0, 0.00359, 0.000424, 0, 0.000424, 0, 0, 0.001333, 0, 0, 0, 0.001326, 0, 0.002336, 0.002408, 0.000424, 0.002902, 0, 0, 0.002361, 0, 0, 0.004338, 0.00087, 0, 0.001356, 0, 0, 0.00087, 0.003039, 0.000424, 0.006266, 0, 0, 0.002405, 0.003591, 0.000426, 0.002357, 0.000435, 0.001358, 0.004681, 0.002691, 0, 0.00902, 0.000866, 0, 0.002355, 0.000871, 0, 0.004251, 0, 0, 0.001805, 0.001847, 0, 0.001323, 0.000867, 0.000867, 0.002418, 0, 0, 0, 0, 0, 0.003326, 0, 0, 0.002368, 0, 0, 0.000423, 0.000424, 0.000424, 0.003706, 0, 0, 0.003546, 0, 0.001893, 0.000919, 0, 0, 0.000426, 0, 0, 0.006524, 0, 0, 0.001955, 0.000423, 0.00087, 0.003471, 0, 0, 0.000428, 0.000435, 0.000435, 0.000881, 0, 0, 0.003216, 0.002908, 0, 0.007752, 0, 0, 0.002305, 0, 0, 0.006781, 0.003127, 0, 0.003127, 0, 0.001349, 0.00305, 0, 0, 0.003765, 0.000428, 0, 0.005602, 0, 0, 0.000866, 0.000868, 0.000868, 0.000866, 0, 0, 0.000426, 0, 0, 0.00435, 0, 0, 0.006003, 0.000871, 0, 0.000428, 0, 0, 0.003064, 0, 0, 0.00088, 0, 0, 0, 0.000423, 0, 0.004115, 0, 0, 0.005536, 0.000426, 0, 0.000423, 0, 0.000423, 0, 0, 0, 0.002876, 0, 0, 0.000425, 0, 0, 0.002879, 0, 0, 0.002351, 0, 0, 0.002352, 0.001942, 0, 0.001942, 0, 0, 0, 0, 0, 0.00133, 0, 0, 0.000877, 0, 0, 0.004291, 0, 0, 0.006154, 0, 0, 0.005701, 0.000424, 0, 0.004187, 0.000423, 0, 0.00088, 0, 0, 0.002352, 0, 0, 0.003438, 0, 0, 0.00684, 0, 0, 0.00937, 0, 0, 0, 0, 0, 0.000916, 0, 0.000428, 0.006328, 0, 0, 0.001443, 0, 0, 0.001935, 0, 0, 0.003471, 0, 0, 0.00235, 0, 0, 0.005219, 0, 0, 0.001851, 0, 0,0.005637, 0, 0, 0, 0, 0, 0.002673, 0, 0.000428, 0, 0, 0, 0.001405, 0, 0, 0.002335, 0, 0, 0, 0, 0, 0.007385, 0, 0, 0, 0, 0, 0.000871, 0,0, 0.003282, 0, 0, 0.003464, 0.000423, 0, 0.001809, 0.000426, 0, 0.001938, 0, 0, 0.012756, 0.000428, 0, 0.002343, 0.000427, 0, 0.004977, 0, 0, 0.001794, 0, 0, 0, 0, 0, 0.001827, 0, 0, 0.002322, 0, 0, 0, 0, 0, 0, 0, 0, 0.002388, 0, 0, 0.003643)
                         
                         
                         if(length(internode)!=5){
                           print("Internode distance list not correct")
                           return
                         }
                         Mu_<- 1/2/(a*Pi_T*Pi_C + b*Pi_T*Pi_A + c*Pi_T*Pi_G +d*Pi_C*Pi_A + e*Pi_C*Pi_G + f*Pi_A*Pi_G)
                         ##Construct Q Matrix
                         Q<-matrix(nrow=4,ncol=4)
                         Q[1,1]<-((-a)*Pi_C) - (b*Pi_A) - (c*Pi_G)
                         Q[1,2]<-a*Pi_C
                         Q[1,3]<-b*Pi_A
                         Q[1,4]<-c*Pi_G
                         Q[2,1]<-a*Pi_T
                         Q[2,2]<-((-a)*Pi_T) - (d*Pi_A) - (e*Pi_G)
                         Q[2,3]<-d*Pi_A
                         Q[2,4]<-e*Pi_G
                         Q[3,1]<-b*Pi_T
                         Q[3,2]<-d*Pi_C
                         Q[3,3]<-((-b)*Pi_T) - (d*Pi_C) - (f*Pi_G)
                         Q[3,4]<-f*Pi_G
                         Q[4,1]<-c*Pi_T
                         Q[4,2]<-e*Pi_C
                         Q[4,3]<-f*Pi_A
                         Q[4,4]<-((-c)*Pi_T) - (e*Pi_C) - (f*Pi_A)
                         Q<-Mu_*Q
                         
                         #Vectorize the base frequencies
                         frequ<-c(Pi_T,Pi_C,Pi_A,Pi_G)
                         
                         #Obtain the eigenvalues and vectors
                         evects<-eigen(Q)
                         evalues<-evects$values
                         evectors<-evects$vectors
                         
                         #Reorder in ascending order, swap rows.
                         evalues<-evalues[c(4,3,2,1)] #same as mathematica
                         evectors<-evectors[c(4,3,2,1),c(4,3,2,1)]
                         
                         
                         #tev<-t(evectors) #depriciated
                         tev<-(evectors)
                         #Get inverse
                         itev<-solve(tev)
                         
                         #Internal function to evaluate lamda
                         evalLambda<-function(lamda){
                           p<-list()
                           p<-array(,dim=c(5,4,4))
                           for(v in 1:length(internode)){
                             p[v,,]<-(tev %*% (diag(exp(evalues*lamda*internode[v]))%*%itev))
                           }
                           
                           correct<-0
                           wrong1<-0
                           wrong2<-0
                           
                           for(original_character in 1:4){
                             for(internode_character in 1:4){
                               for(leaf_character_1 in 1:4){
                                 for(leaf_character_2 in 1:4){
                                   if(leaf_character_1!=leaf_character_2){
                                     correct <- correct+(frequ[original_character]*
                                                           p[5,original_character, internode_character]*
                                                           p[1,original_character, leaf_character_1]*
                                                           p[2,original_character, leaf_character_1]*
                                                           p[3,internode_character, leaf_character_2]*
                                                           p[4,internode_character, leaf_character_2])
                                     wrong1 <- wrong1+(frequ[original_character]*
                                                         p[5,original_character, internode_character]*
                                                         p[1,original_character, leaf_character_1]*
                                                         p[2,original_character, leaf_character_2]*
                                                         p[3,internode_character, leaf_character_1]*
                                                         p[4,internode_character, leaf_character_2])
                                     wrong2<-wrong2+ (frequ[original_character]*
                                                        p[5,original_character, internode_character]*
                                                        p[1,original_character, leaf_character_1]*
                                                        p[2,original_character, leaf_character_2]*
                                                        p[3,internode_character, leaf_character_2]*
                                                        p[4,internode_character, leaf_character_1])
                                   }
                                 }
                               }
                             }
                           }
                           all<-c(correct,wrong1,wrong2)
                           return(all)
                         }
                         #Initialize blanks
                         eYsum <- 0
                         eX1sum <- 0
                         eX2sum <- 0
                         eY2sum <- 0
                         eX12sum <- 0
                         eX22sum <- 0
                         eX1Ysum <- 0
                         eX2Ysum <- 0
                         eX1X2sum <- 0
                         
                         for(lmbda in rate_vector){
                           all<-evalLambda(lmbda)
                           y<-all[1]
                           x1<-all[2]
                           x2<-all[3]	
                           eYsum<-eYsum+y
                           eX1sum<-eX1sum+x1
                           eX2sum<-eX2sum+x2
                           
                           eY2sum<-eY2sum+(y^2)
                           eX12sum<-eX12sum+(x1^2)
                           eX22sum<-eX22sum+(x2^2)
                           
                           eX1Ysum<-eX1Ysum+(x1*y)
                           eX2Ysum<-eX2Ysum+(x2*y)
                           eX1X2sum<-eX1X2sum+(x1*x2)
                         }
                         
                         Mu_1 <- eYsum - eX1sum
                         Mu_2 <- eYsum - eX2sum
                         
                         
                         Sigma_1 <- sqrt(eX1sum + eYsum - eX12sum - eY2sum + 2*eX1Ysum)
                         Sigma_2 <-sqrt(eX2sum + eYsum - eX22sum - eY2sum + 2*eX2Ysum)
                         Rho_<- (-eX1X2sum + eX1Ysum + eX2Ysum + eYsum - eY2sum)/(Sigma_1*Sigma_2)
                         
                         #Internal function for integration
                         FofT<-function(t){
                           F1ofT=((1 / Sigma_1) * dnorm((t - Mu_1)/ Sigma_1)*pnorm(Rho_*(t - Mu_1)/(Sigma_1* sqrt(1 - Rho_*Rho_)) - (t - Mu_2)/(Sigma_2* sqrt(1 - Rho_*Rho_))))
                           F2ofT=((1 / Sigma_2) * dnorm((t - Mu_2)/ Sigma_2)*pnorm(Rho_*(t - Mu_2)/(Sigma_2* sqrt(1 - Rho_*Rho_)) - (t - Mu_1)/(Sigma_1* sqrt(1 - Rho_*Rho_))))
                           return(F1ofT+F2ofT)
                         }
                         
                         princtree<-integrate(FofT, -Inf, -.5)
                         prpolytomy = integrate(FofT, -.5, .5)
                         prcortree  = integrate(FofT, .5, Inf)
                         
                         print(paste0("Probablility Correct: ",prcortree$value))
                         print(paste0("Probability Incorrect: ",princtree$value))
                         print(paste0("Probability Polytomy: ",prpolytomy$value))
                         return(c(princtree$value,prpolytomy$value,prcortree$value))
}


get.tree<-function(quart,tree){
  as.matrix(tree$tip.label)->drop
  drop[which(!drop[,1]%in%quart),]->prune
  drop.tip(tree,prune)->four.taxa	
  return(four.taxa)
}
bayes.signal.prep<-function(quart,tree){
  
  get.tree(quart,tree)->four.taxa
  combn(quart,2)->get
  knowledge<-matrix()
  for (i in 1:6){
    is.monophyletic(four.taxa,get[,i])->knowledge[i]
  }
  
  length(which(knowledge[1:6]=="TRUE"))->pec.or.quart
  
  if (pec.or.quart==2)	{
    
    ##the following line from Liam Revells phytools blog
    ee<-setNames(four.taxa$edge.length[sapply(1:4,function(x,y) which(y==x), y=four.taxa$edge[,2])],four.taxa$tip.label)
    max(branching.times(four.taxa))-max(ee)->internode
    "internode"->names(internode)
    ##arrange
    which(names(ee)==quart[1])->first
    which(names(ee)==quart[2])->second
    which(names(ee)==quart[3])->third
    which(names(ee)==quart[4])->fourth
    c(ee[first],ee[second],ee[third],ee[fourth], internode)->vector
  }	
  else if (pec.or.quart==1){
    ee<-setNames(four.taxa$edge.length[sapply(1:4,function(x,y) which(y==x), y=four.taxa$edge[,2])],four.taxa$tip.label)
    ##arrange to get to what the internode is and where to add the BL to T1.
    rev(sort(ee))->ee2
    ee2[1]-ee2[2]->internode
    "internode"->names(internode)
    ee2+c(internode,0,0,0)->newee
    which(names(newee)==quart[1])->first
    which(names(newee)==quart[2])->second
    which(names(newee)==quart[3])->third
    which(names(newee)==quart[4])->fourth
    c(newee[first], newee[second], newee[third], newee[fourth], internode)->vector	
  }	
  
  return(vector)
}



##these are the same inputs as the allmodel.signal.noise, users will use this and save the output for plotting
post.su<-function(a,b,c,d,e,f,Pi_T,Pi_C,Pi_A,Pi_G, rate_vector,quart,tree)
{
  
  ###first get your internodes
  matrix(ncol=5)->stored_ints
  for (i in 1:length(tree))
  {
    bayes.signal.prep(quart,tree[[i]])-> temp
    rbind(stored_ints,temp)->stored_ints
    
  }
  stored_ints[2:length(stored_ints[,1]),]->stored_ints
  length(stored_ints[,1])->loop.length
  matrix(ncol=length(stored_ints[,1]),nrow=3)-> quart.probs
  #foreach(i=1:loop.length, .combine=cbind)%dopar%
  for (i in 2:length(stored_ints[,1]))
  {
    allmodel.signal.noise (a,b,c,d,e,f, stored_ints[i,],Pi_T,Pi_C,Pi_A,Pi_G, rate_vector)-> temp2 #quart.probs[,i]
    #rbind(quart.probs,temp2)->quart.probs#[i,]
    temp2->quart.probs[,i]
  }
  t(quart.probs)->qp2
  cbind(qp2,stored_ints)->final
  return(qp2)	
  
}

##### User function. foreach does not work with downstream manipulations of objects well, so this takes the output of the core post.su function and adds the internode lengths back to have one nice result object
#' @export
su.bayes<-function(a,b,c,d,e,f,Pi_T,Pi_C,Pi_A,Pi_G, rate_vector,quart,tree){
  post.su(a,b,c,d,e,f,Pi_T,Pi_C,Pi_A,Pi_G, rate_vector,quart,tree)->final
  t(final)->qp2
  matrix(ncol=5)->stored_ints
  for (i in 1:length(tree))
  {
    bayes.signal.prep(quart,tree[[i]])-> temp
    rbind(stored_ints,temp)->stored_ints
  }
  cbind(t(qp2),stored_ints[2:length(stored_ints[,1]),])->final.result
  return(final.result)	
  
}


###This will either plot the Quartet internode probs with their internode, or else the violin plots o look at density another way
#' @export
plotPosterior<-function(final, plotType="QIPs")
{
  as.data.frame(final)->final2
  ##Experimental
  final2<-final2[2:nrow(final2),]
  ##  
  dim(final)->ll
  ll[1]->up
  final2[2:up,]->final22
  x    <- as.numeric(as.character(final22[,8]))
  y1    <- as.numeric(as.character(final22[,3]))
  y2    <- as.numeric(as.character(final22[,2])) #polytomy
  y3    <- as.numeric(as.character(final22[,1]))
  
  if (plotType=="QIPs")	{
    p1<-ggplot(final22,aes(x=x,y=y1)) + stat_binhex(colour="white",na.rm=TRUE)+ xlab("internode length") + ylab("QIRP") + scale_fill_gradientn(colours=c("green1","red"),name = "Frequency",na.value=NA)+ theme_bw()
    p2<-ggplot(final22,aes(x=x,y=y2)) + stat_binhex(colour="white",na.rm=TRUE)+ xlab("internode length") + ylab("QIPP") + scale_fill_gradientn(colours=c("green1","red"),name = "Frequency",na.value=NA)+ theme_bw()
    p3<-ggplot(final22,aes(x=x,y=y3)) + stat_binhex(colour="white",na.rm=TRUE)+ xlab("internode length") + ylab("QIHP") + scale_fill_gradientn(colours=c("green1","red"),name = "Frequency",na.value=NA)+ theme_bw()
    grid.arrange(p1, p2, p3, ncol=1, nrow =3)
  } else if (plotType=="violin"){
    c(y1,y2,y3)->stacks
    length(y1)->set
    rep("QIRP",set)->Qirp
    rep("QIPP",set)->Qipp
    rep("QIHP",set)->Qihp
    c(Qirp,Qipp,Qihp)->c2
    cbind(stacks,c2)->newy
    colnames(newy)<-c("Probability","Analysis")
    rep(x,3)->internodes
    #colnames(internodes)<-"internode"
    cbind(internodes,newy)->data
    as.data.frame(data)->data
    as.numeric(as.character(data[,1]))->data[,1]
    as.numeric(as.character(data[,2]))->data[,2]
    Analysis<-data[,"Analysis"]
    Probability<-data[,"Probability"]
    p<-ggplot(data, aes(x= Analysis, y= Probability, fill=Analysis)) 
    p + geom_violin(trim=FALSE)+scale_fill_manual(values=c("firebrick","deepskyblue3","seagreen")) + geom_boxplot(width=0.1, fill= "aliceblue")		
    
  }}






###############################################################################
## Modernization overrides (performance + tidyverse binding + interactive HTML)
## - Drop-in overrides: later definitions replace earlier ones in R.
## - Goal: speed up execution without changing the underlying algorithm.
## - Compatibility: preserve key object shapes and RNG side-effects where relevant.
##
## Usage:
##   - Set RUN_SELF_TEST <- TRUE once to validate “before vs after” consistency.
##   - After PASS, set it back to FALSE for normal runs.
###############################################################################

RUN_SELF_TEST <- FALSE  # Set TRUE once to run internal consistency checks.

## ---------------------------------------------------------------------------
## (0) Baseline capture (computed BEFORE overriding functions)
##     This allows in-file regression checks without creating a second script.
## ---------------------------------------------------------------------------
if (RUN_SELF_TEST) {
  set.seed(123)
  rv_test <- rep(c(0.1, 0.2, 0.2, 0.3), 200)
  t_test  <- 1.0
  t0_test <- 0.2
  s_test  <- 4
  tree_test <- ape::rtree(6)
  
  base_site <- site.summer(rv_test, t_test)
  base_ipg2_dim <- dim(inform.profile.generator2(rv_test, tree_test))
  
  base_App <- Approximator(t_test, t0_test, rv_test, s_test)
  
  set.seed(999)
  base_Ex  <- ExMaker(t_test, t0_test, rv_test, s_test)
  
  internode_test <- c(0.1, 0.1, 0.1, 0.1, 0.05)
  ratevec_test <- c(0.01, 0.02, 0.02, 0.03)
  base_all <- allmodel.signal.noise(
    1,1,1,1,1,1, internode_test,
    0.25,0.25,0.25,0.25, ratevec_test
  )
}

## ---------------------------------------------------------------------------
## (1) Exact site-pattern folding helper (NO rounding / NO approximation)
##     This is mathematically exact whenever the downstream logic is summation-
##     based over sites (which it is for the functions we optimize here).
## ---------------------------------------------------------------------------
.fold_rates_exact <- function(rateVector) {
  rv <- as.numeric(rateVector)
  rv <- rv[!is.na(rv)]
  u <- unique(rv)
  idx <- match(rv, u)
  w <- tabulate(idx, nbins = length(u))
  list(lambda = u, weight = w, n = length(rv))
}

## ---------------------------------------------------------------------------
## (2) Performance overrides (algorithm unchanged)
## ---------------------------------------------------------------------------

## (2a) site.summer: strict vectorization of the Townsend 2007 expression
site.summer <- function(rate.vector, time) {
  rv <- as.numeric(rate.vector)
  sum(16 * rv * rv * time * exp(-4 * rv * time))
}

## (2b) inform.profile.generator2: preserve original 1×k matrix output shape
inform.profile.generator2 <- function(use.rates, tree) {
  btimes2 <- sort(c(0, branching.times(tree)))
  out <- vapply(btimes2, function(bt) site.summer(use.rates, bt), numeric(1))
  matrix(out, nrow = 1)  # IMPORTANT: keep 1×k matrix for downstream compatibility
}

## (2c) inform.profile.generator: keep original return structure (2×k matrix)
inform.profile.generator <- function(rate.vector, tree) {
  btimes2 <- sort(c(0, branching.times(tree)))
  out <- vapply(btimes2, function(bt) site.summer(rate.vector, bt), numeric(1))
  rbind(btimes2, out)
}

## (2d) Approximator: exact folding + vectorized sums; keep RNG consumption
Approximator <- function(t, t0, rateVector, s) {
  fr <- .fold_rates_exact(rateVector)
  lam <- fr$lambda
  w   <- fr$weight
  n   <- fr$n
  
  npnl  <- pnl(lam, t, t0, s)
  npsnr <- psnr(lam, t, t0, s)
  
  Ey    <- sum(w * npsnr)
  Ex1   <- sum(w * npnl)
  Ex2   <- sum(w * (npnl * npnl))
  Ey2   <- sum(w * (npsnr * npsnr))
  EX1Y  <- sum(w * (npsnr * npnl))
  ESQ   <- sum(w * (npsnr * sqrt(npnl)))  # Retains original NaN behavior if npnl < 0
  
  Ex <- Ex1 + sqrt(Ex1 / pi)
  Expectation <- Ey - Ex
  
  variancey <- Ey - Ey2
  variancex <- ((pi - 1) / pi) * Ex1 - Ex2
  variance  <- variancey + variancex - 2 * EX1Y - (2 / sqrt(pi)) * ESQ
  
  ## Preserve RNG side-effect: original code draws rnorm(n, ...) even though unused.
  rnorm(n, mean = Expectation, sd = sqrt(variance))
  
  princtree  <- pnorm(-0.5, mean = Expectation, sd = sqrt(variance))
  prpolytomy <- pnorm(0.5, mean = Expectation, sd = sqrt(variance)) - princtree
  prcortree  <- 1 - pnorm(0.5, mean = Expectation, sd = sqrt(variance))
  
  values <- c(prcortree, prpolytomy, princtree)
  names(values) <- c("Probabilty Correct", "Probability Polytomy", "Probability Incorrect")
  values
}

## (2e) Approximator.lite: same as above but returns only P(correct)
Approximator.lite <- function(t, t0, rateVector, s) {
  fr <- .fold_rates_exact(rateVector)
  lam <- fr$lambda
  w   <- fr$weight
  n   <- fr$n
  
  npnl  <- pnl(lam, t, t0, s)
  npsnr <- psnr(lam, t, t0, s)
  
  Ey    <- sum(w * npsnr)
  Ex1   <- sum(w * npnl)
  Ex2   <- sum(w * (npnl * npnl))
  Ey2   <- sum(w * (npsnr * npsnr))
  EX1Y  <- sum(w * (npsnr * npnl))
  ESQ   <- sum(w * (npsnr * sqrt(npnl)))
  
  Ex <- Ex1 + sqrt(Ex1 / pi)
  Expectation <- Ey - Ex
  
  variancey <- Ey - Ey2
  variancex <- ((pi - 1) / pi) * Ex1 - Ex2
  variance  <- variancey + variancex - 2 * EX1Y - (2 / sqrt(pi)) * ESQ
  
  ## Preserve RNG side-effect.
  rnorm(n, mean = Expectation, sd = sqrt(variance))
  
  1 - pnorm(0.5, mean = Expectation, sd = sqrt(variance))
}

## (2f) ExMaker: vectorized runif with identical RNG stream (vs n×runif(1))
ExMaker <- function(t, t0, rateVector, s) {
  rv <- as.numeric(rateVector)
  rv <- rv[!is.na(rv)]
  n  <- length(rv)
  
  npro  <- prother(rv, t, t0, s)
  npsnr <- psnr(rv, t, t0, s)
  npnl  <- pnl(rv,  t, t0, s)
  
  u <- runif(n, min = 0, max = 1)
  
  ## Replicates the original if/else chain outcome.
  snr   <- (u >= npro) & (u < (npro + npsnr))
  nwone <- (u >= (npro + npsnr)) & (u < (npro + npsnr + npnl))
  
  snr_sum   <- sum(snr)
  nwone_sum <- sum(nwone)
  
  ## Preserve the effective behavior of the current implementation:
  ## do NOT “fix” potential logic issues here, to avoid changing results.
  nwtwo_sum <- 0
  
  wronger <- max(c(nwone_sum, nwtwo_sum))
  n + 1 + snr_sum - wronger
}

## ---------------------------------------------------------------------------
## (3) allmodel.signal.noise: algebraic speedup + exact folding
##     - Replaces 4-level leaf loops with an equivalent closed-form reduction:
##       sum_{l1!=l2} A[l1]*B[l2] = sum(A)*sum(B) - sum(A*B)
##     - Adds exact folding over rate_vector to reduce repeated evalLambda calls.
## ---------------------------------------------------------------------------
allmodel.signal.noise <- function(a,b,c,d,e,f,internode,Pi_T,Pi_C,Pi_A,Pi_G, rate_vector) {
  
  if (length(internode) != 5) {
    print("Internode distance list not correct")
    return()
  }
  
  Mu_<- 1/2/(a*Pi_T*Pi_C + b*Pi_T*Pi_A + c*Pi_T*Pi_G +
               d*Pi_C*Pi_A + e*Pi_C*Pi_G + f*Pi_A*Pi_G)
  
  ## Construct Q matrix
  Q <- matrix(nrow=4, ncol=4)
  Q[1,1] <- (-a*Pi_C) - (b*Pi_A) - (c*Pi_G)
  Q[1,2] <-  a*Pi_C
  Q[1,3] <-  b*Pi_A
  Q[1,4] <-  c*Pi_G
  Q[2,1] <-  a*Pi_T
  Q[2,2] <- (-a*Pi_T) - (d*Pi_A) - (e*Pi_G)
  Q[2,3] <-  d*Pi_A
  Q[2,4] <-  e*Pi_G
  Q[3,1] <-  b*Pi_T
  Q[3,2] <-  d*Pi_C
  Q[3,3] <- (-b*Pi_T) - (d*Pi_C) - (f*Pi_G)
  Q[3,4] <-  f*Pi_G
  Q[4,1] <-  c*Pi_T
  Q[4,2] <-  e*Pi_C
  Q[4,3] <-  f*Pi_A
  Q[4,4] <- (-c*Pi_T) - (e*Pi_C) - (f*Pi_A)
  Q <- Mu_ * Q
  
  frequ <- c(Pi_T,Pi_C,Pi_A,Pi_G)
  
  ## Eigen decomposition
  evects  <- eigen(Q)
  evalues <- evects$values
  evectors<- evects$vectors
  
  ## Match original ordering
  evalues  <- evalues[c(4,3,2,1)]
  evectors <- evectors[c(4,3,2,1), c(4,3,2,1)]
  
  tev  <- evectors
  itev <- solve(tev)
  
  evalLambda <- function(lamda) {
    p <- array(0, dim = c(5,4,4))
    for (v in 1:length(internode)) {
      p[v,,] <- (tev %*% (diag(exp(evalues * lamda * internode[v])) %*% itev))
    }
    
    correct <- 0
    wrong1  <- 0
    wrong2  <- 0
    
    for (oc in 1:4) {
      for (ic in 1:4) {
        
        ## correct
        A <- p[1, oc, ] * p[2, oc, ]
        B <- p[3, ic, ] * p[4, ic, ]
        corr_term <- sum(A) * sum(B) - sum(A * B)
        
        ## wrong1
        C <- p[1, oc, ] * p[3, ic, ]
        D <- p[2, oc, ] * p[4, ic, ]
        w1_term <- sum(C) * sum(D) - sum(C * D)
        
        ## wrong2
        E <- p[1, oc, ] * p[4, ic, ]
        F <- p[2, oc, ] * p[3, ic, ]
        w2_term <- sum(E) * sum(F) - sum(E * F)
        
        wt <- frequ[oc] * p[5, oc, ic]
        correct <- correct + wt * corr_term
        wrong1  <- wrong1  + wt * w1_term
        wrong2  <- wrong2  + wt * w2_term
      }
    }
    
    c(correct, wrong1, wrong2)
  }
  
  ## Accumulators
  eYsum <- 0; eX1sum <- 0; eX2sum <- 0
  eY2sum <- 0; eX12sum <- 0; eX22sum <- 0
  eX1Ysum <- 0; eX2Ysum <- 0; eX1X2sum <- 0
  
  ## Exact folding over rate vector
  fr <- .fold_rates_exact(rate_vector)
  lam_u <- fr$lambda
  w_u   <- fr$weight
  
  for (k in seq_along(lam_u)) {
    lmbda <- lam_u[k]
    wgt   <- w_u[k]
    
    all <- evalLambda(lmbda)
    y  <- all[1]
    x1 <- all[2]
    x2 <- all[3]
    
    eYsum  <- eYsum  + wgt * y
    eX1sum <- eX1sum + wgt * x1
    eX2sum <- eX2sum + wgt * x2
    
    eY2sum  <- eY2sum  + wgt * (y^2)
    eX12sum <- eX12sum + wgt * (x1^2)
    eX22sum <- eX22sum + wgt * (x2^2)
    
    eX1Ysum  <- eX1Ysum  + wgt * (x1*y)
    eX2Ysum  <- eX2Ysum  + wgt * (x2*y)
    eX1X2sum <- eX1X2sum + wgt * (x1*x2)
  }
  
  Mu_1 <- eYsum - eX1sum
  Mu_2 <- eYsum - eX2sum
  
  Sigma_1 <- sqrt(eX1sum + eYsum - eX12sum - eY2sum + 2*eX1Ysum)
  Sigma_2 <- sqrt(eX2sum + eYsum - eX22sum - eY2sum + 2*eX2Ysum)
  Rho_<- (-eX1X2sum + eX1Ysum + eX2Ysum + eYsum - eY2sum)/(Sigma_1*Sigma_2)
  
  FofT <- function(t) {
    F1 <- (1 / Sigma_1) * dnorm((t - Mu_1)/Sigma_1) *
      pnorm(Rho_*(t - Mu_1)/(Sigma_1*sqrt(1 - Rho_*Rho_)) -
              (t - Mu_2)/(Sigma_2*sqrt(1 - Rho_*Rho_)))
    F2 <- (1 / Sigma_2) * dnorm((t - Mu_2)/Sigma_2) *
      pnorm(Rho_*(t - Mu_2)/(Sigma_2*sqrt(1 - Rho_*Rho_)) -
              (t - Mu_1)/(Sigma_1*sqrt(1 - Rho_*Rho_)))
    F1 + F2
  }
  
  princtree  <- integrate(FofT, -Inf, -.5)
  prpolytomy <- integrate(FofT, -.5,  .5)
  prcortree  <- integrate(FofT,  .5,  Inf)
  
  print(paste0("Probablility Correct: ", prcortree$value))
  print(paste0("Probability Incorrect: ", princtree$value))
  print(paste0("Probability Polytomy: ", prpolytomy$value))
  c(princtree$value, prpolytomy$value, prcortree$value)
}

## ---------------------------------------------------------------------------
## (4) tidyverse binding overrides: remove rbind-in-loop in post.su / su.bayes
##     - Preserve legacy behavior: keep the first probability row as NA
##       (your downstream plotting code skips the first row anyway).
## ---------------------------------------------------------------------------
post.su <- function(a,b,c,d,e,f,Pi_T,Pi_C,Pi_A,Pi_G, rate_vector, quart, tree) {
  
  stored_ints <- purrr::map(tree, ~ bayes.signal.prep(quart, .x)) %>%
    purrr::map(~ matrix(.x, nrow = 1)) %>%
    dplyr::bind_rows() %>%
    as.matrix()
  
  loop.length <- nrow(stored_ints)
  
  quart.probs <- matrix(NA_real_, nrow = 3, ncol = loop.length)
  for (i in 2:loop.length) {
    quart.probs[, i] <- allmodel.signal.noise(
      a,b,c,d,e,f, stored_ints[i,], Pi_T,Pi_C,Pi_A,Pi_G, rate_vector
    )
  }
  
  t(quart.probs)
}

su.bayes <- function(a,b,c,d,e,f,Pi_T,Pi_C,Pi_A,Pi_G, rate_vector, quart, tree) {
  qp2 <- post.su(a,b,c,d,e,f,Pi_T,Pi_C,Pi_A,Pi_G, rate_vector, quart, tree)
  
  stored_ints <- purrr::map(tree, ~ bayes.signal.prep(quart, .x)) %>%
    purrr::map(~ matrix(.x, nrow = 1)) %>%
    dplyr::bind_rows() %>%
    as.matrix()
  
  cbind(qp2, stored_ints)
}

## ---------------------------------------------------------------------------
## (5) Interactive HTML visualization (new functions; legacy ggplot remains)
##     - Hover tooltips for detailed values.
##     - Legend click-to-toggle (built-in plotly behavior).
## ---------------------------------------------------------------------------
plotPosterior_html <- function(final, plotType = c("QIPs","violin"),
                               file = NULL, selfcontained = TRUE) {
  
  plotType <- match.arg(plotType)
  final2 <- as.data.frame(final)
  
  ## Preserve your legacy skip-first-row behavior
  final2 <- final2[2:nrow(final2), , drop = FALSE]
  up <- nrow(final2)
  dat <- final2[2:up, , drop = FALSE]
  
  x  <- as.numeric(as.character(dat[,8]))
  y1 <- as.numeric(as.character(dat[,3]))  # QIRP
  y2 <- as.numeric(as.character(dat[,2]))  # QIPP
  y3 <- as.numeric(as.character(dat[,1]))  # QIHP
  
  if (plotType == "QIPs") {
    p1 <- plot_ly(x = x, y = y1, type = "histogram2d",
                  hovertemplate = "internode=%{x}<br>QIRP=%{y}<br>count=%{z}<extra></extra>") %>%
      layout(xaxis = list(title = "internode length"),
             yaxis = list(title = "QIRP"))
    p2 <- plot_ly(x = x, y = y2, type = "histogram2d",
                  hovertemplate = "internode=%{x}<br>QIPP=%{y}<br>count=%{z}<extra></extra>") %>%
      layout(xaxis = list(title = "internode length"),
             yaxis = list(title = "QIPP"))
    p3 <- plot_ly(x = x, y = y3, type = "histogram2d",
                  hovertemplate = "internode=%{x}<br>QIHP=%{y}<br>count=%{z}<extra></extra>") %>%
      layout(xaxis = list(title = "internode length"),
             yaxis = list(title = "QIHP"))
    
    w <- subplot(p1, p2, p3, nrows = 3, shareX = TRUE, titleY = TRUE) %>%
      layout(showlegend = FALSE)
    
  } else {
    df <- tibble(
      internode = rep(x, 3),
      Probability = c(y1, y2, y3),
      Analysis = factor(rep(c("QIRP","QIPP","QIHP"), each = length(x)),
                        levels = c("QIRP","QIPP","QIHP"))
    )
    
    w <- plot_ly(df, x = ~Analysis, y = ~Probability, type = "violin",
                 color = ~Analysis, box = list(visible = TRUE),
                 points = "outliers",
                 hovertemplate = "Type=%{x}<br>Prob=%{y}<extra></extra>") %>%
      layout(xaxis = list(title = ""),
             yaxis = list(title = "Probability"))
  }
  
  if (!is.null(file)) {
    htmlwidgets::saveWidget(w, file = file, selfcontained = selfcontained)
  }
  w
}

graph.signal.noise_html <- function(currentprobdist, rateVector,
                                    file = NULL, selfcontained = TRUE) {
  n <- length(rateVector)
  max.bound <- 2*n + 1
  start.right <- n + 2
  poly <- n + 1
  
  m <- as.numeric(currentprobdist)
  polytomy <- m[poly]
  blue.side <- m[start.right:max.bound]
  grey.side <- m[1:n]
  
  plotter <- c(grey.side, polytomy, blue.side)
  x <- seq_along(plotter)
  group <- c(rep("wrong", length(grey.side)), "polytomy", rep("correct", length(blue.side)))
  
  df <- tibble(x = x, freq = plotter, group = group)
  
  w <- plot_ly(df, x = ~x, y = ~freq, type = "bar", color = ~group,
               customdata = ~group,
               hovertemplate = "class=%{customdata}<br>x=%{x}<br>freq=%{y}<extra></extra>") %>%
    layout(xaxis = list(title = "Signal/Noise index"),
           yaxis = list(title = "Frequency"))
  
  if (!is.null(file)) {
    htmlwidgets::saveWidget(w, file = file, selfcontained = selfcontained)
  }
  w
}

## ---------------------------------------------------------------------------
## (6) Post-override self-test (compare to captured baselines)
## ---------------------------------------------------------------------------
if (RUN_SELF_TEST) {
  set.seed(123)
  rv_test <- rep(c(0.1, 0.2, 0.2, 0.3), 200)
  t_test  <- 1.0
  t0_test <- 0.2
  s_test  <- 4
  tree_test <- ape::rtree(6)
  
  stopifnot(isTRUE(all.equal(base_site, site.summer(rv_test, t_test), tolerance = 1e-12)))
  stopifnot(identical(base_ipg2_dim, dim(inform.profile.generator2(rv_test, tree_test))))
  stopifnot(isTRUE(all.equal(base_App, Approximator(t_test, t0_test, rv_test, s_test), tolerance = 1e-12)))
  
  set.seed(999)
  stopifnot(identical(base_Ex, ExMaker(t_test, t0_test, rv_test, s_test)))
  
  internode_test <- c(0.1, 0.1, 0.1, 0.1, 0.05)
  ratevec_test <- c(0.01, 0.02, 0.02, 0.03)
  new_all <- allmodel.signal.noise(1,1,1,1,1,1, internode_test, 0.25,0.25,0.25,0.25, ratevec_test)
  stopifnot(isTRUE(all.equal(base_all, new_all, tolerance = 1e-12)))
  
  message("PASS: modernization overrides match baseline within tolerance.")
}




informativeness.profile_html <- function(rate.vector, tree, codon = "FALSE",
                                         file = NULL, selfcontained = TRUE) {
  btimes <- sort(c(0, branching.times(tree)))
  
  if (codon == "FALSE") {
    y <- vapply(btimes, function(bt) site.summer(rate.vector, bt), numeric(1))
    df <- data.frame(time = btimes, PI = y)
    
    w <- plotly::plot_ly(
      df, x = ~time, y = ~PI, type = "scatter", mode = "lines",
      hovertemplate = "time=%{x}<br>PI=%{y}<extra></extra>"
    ) %>%
      plotly::layout(xaxis = list(title = "Time from present"),
                     yaxis = list(title = "Phylogenetic informativeness"))
    
  } else {
    pos1 <- rate.vector[seq(1, length(rate.vector), 3)]
    pos2 <- rate.vector[seq(2, length(rate.vector), 3)]
    pos3 <- rate.vector[seq(3, length(rate.vector), 3)]
    
    y1 <- vapply(btimes, function(bt) site.summer(pos1, bt), numeric(1))
    y2 <- vapply(btimes, function(bt) site.summer(pos2, bt), numeric(1))
    y3 <- vapply(btimes, function(bt) site.summer(pos3, bt), numeric(1))
    
    df <- data.frame(time = btimes, pos1 = y1, pos2 = y2, pos3 = y3)
    
    w <- plotly::plot_ly(df, x = ~time, y = ~pos1, type = "scatter", mode = "lines",
                         name = "pos1",
                         hovertemplate = "time=%{x}<br>pos1=%{y}<extra></extra>") %>%
      plotly::add_lines(y = ~pos2, name = "pos2",
                        hovertemplate = "time=%{x}<br>pos2=%{y}<extra></extra>") %>%
      plotly::add_lines(y = ~pos3, name = "pos3",
                        hovertemplate = "time=%{x}<br>pos3=%{y}<extra></extra>") %>%
      plotly::layout(xaxis = list(title = "Time from present"),
                     yaxis = list(title = "Phylogenetic informativeness"))
  }
  
  if (!is.null(file)) {
    htmlwidgets::saveWidget(w, file = file, selfcontained = selfcontained)
  }
  w
}

## ---------------------------------------------------------------------------
## Multi-locus interactive PI profile (HTML)
## - One trace per locus (gene/partition).
## - Hover shows locus name + peak metrics.
## - Controls: sort strategy + time-of-interest + show top N.
## - UX upgrades:
##   * Sticky control panel (always visible)
##   * Page scroll enabled
##   * Plotly: disable scroll wheel zoom (avoid accidental zoom)
##   * Plotly: default drag = pan (less mis-click)
## - Algorithm unchanged: PI(t) = sum_i 16*lambda_i^2*t*exp(-4*lambda_i*t)
## ---------------------------------------------------------------------------

## Fast exact PI(t) with rate folding (unique values + multiplicities).
## This is mathematically identical to summing over all sites.
.pi_at_time_folded <- function(u_rates, u_counts, t) {
  if (t <= 0) return(0)
  16 * t * sum(u_counts * (u_rates^2) * exp(-4 * u_rates * t))
}

## Compute PI curve for one locus on a common time grid.
.pi_curve_one_locus <- function(rate_vec, times) {
  rate_vec <- as.numeric(rate_vec)
  tab <- table(rate_vec)
  u_rates  <- as.numeric(names(tab))
  u_counts <- as.numeric(tab)
  vapply(times, function(tt) .pi_at_time_folded(u_rates, u_counts, tt), numeric(1))
}

## Main: multi-locus interactive HTML
## rates_list: named list of numeric vectors (each locus)
informativeness.profile_multi_html <- function(rates_list, tree,
                                               times = NULL,
                                               file = NULL,
                                               selfcontained = TRUE,
                                               default_top_n = 50) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Please install.packages('plotly')")
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) stop("Please install.packages('htmlwidgets')")
  if (!requireNamespace("htmltools", quietly = TRUE)) stop("Please install.packages('htmltools')")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Please install.packages('jsonlite')")
  
  if (!is.list(rates_list) || length(rates_list) < 2) {
    stop("rates_list must be a named list with >=2 loci (genes/partitions).")
  }
  if (is.null(names(rates_list)) || any(names(rates_list) == "")) {
    names(rates_list) <- paste0("Locus_", seq_along(rates_list))
  }
  
  ## Default time grid: branching times (matches legacy behavior).
  if (is.null(times)) {
    times <- sort(c(0, branching.times(tree)))
  } else {
    times <- sort(unique(as.numeric(times)))
  }
  
  ## Compute PI curves (exact; algorithm unchanged)
  loci <- names(rates_list)
  curves <- lapply(rates_list, .pi_curve_one_locus, times = times)
  
  ## Peak metrics per locus
  peak_height <- vapply(curves, max, numeric(1), na.rm = TRUE)
  peak_idx    <- vapply(curves, which.max, integer(1))
  peak_time   <- times[peak_idx]
  
  meta <- data.frame(
    locus = loci,
    n_sites = vapply(rates_list, length, integer(1)),
    peak_height = peak_height,
    peak_time = peak_time,
    stringsAsFactors = FALSE
  )
  
  ## Build plotly with one trace per locus (default show top N by peak_height)
  o_default <- order(meta$peak_height, decreasing = TRUE)
  show_n <- min(default_top_n, length(loci))
  visible0 <- rep(FALSE, length(loci))
  visible0[o_default[seq_len(show_n)]] <- TRUE
  
  p <- plotly::plot_ly()
  for (i in seq_along(loci)) {
    nm <- loci[i]
    y  <- curves[[i]]
    p <- plotly::add_lines(
      p, x = times, y = y, name = nm, visible = visible0[i],
      hovertemplate = paste0(
        "locus=", nm,
        "<br>time=%{x}",
        "<br>PI=%{y}",
        "<br>peak_time=", signif(meta$peak_time[i], 6),
        "<br>peak_height=", signif(meta$peak_height[i], 6),
        "<extra></extra>"
      )
    )
  }
  
  ## Add a reference vertical line (default at median time)
  t_ref0 <- stats::median(times)
  
  ## Plot layout + UX defaults
  p <- plotly::layout(
    p,
    height = 650,                 ## fixed height => better page scrolling behavior
    dragmode = "pan",             ## default drag = pan (less accidental zoom)
    xaxis = list(title = "Time from present"),
    yaxis = list(title = "Phylogenetic informativeness"),
    shapes = list(list(
      type = "line", x0 = t_ref0, x1 = t_ref0,
      y0 = 0, y1 = 1, xref = "x", yref = "paper",
      line = list(dash = "dash", width = 2),
      editable = TRUE          ## allow dragging this line
    )),
    edits = list(shapePosition = TRUE),  ## enable dragging shapes
    
    legend = list(orientation = "v")
  )
  
  ## Plotly widget config (safe interactions + export defaults)
  p <- plotly::config(
    p,
    scrollZoom = FALSE,           ## disable mouse-wheel / trackpad zoom
    doubleClick = "reset",        ## double-click resets view
    displaylogo = FALSE,
    editable = TRUE,              ## REQUIRED for dragging the vertical dashed line
    toImageButtonOptions = list(  ## modebar export defaults
      format = "svg",
      filename = "PI_multi_locus",
      scale = 2
    )
  )
  
  ## HTML controls + JS: re-rank and show top N without re-computing curves
  ## Sort modes:
  ## 1) peak_height (desc)
  ## 2) peak_time (desc; "latest peak")
  ## 3) closest_to_ref (asc by |peak_time - t_ref|, then peak_height desc)
  js <- sprintf("
function(el, x){
  var gd = document.getElementById(el.id);
  var meta = %s;

  function order_indices(mode, t_ref){
    var idx = meta.map((d,i)=>i);
    if(mode==='peak_height'){
      idx.sort((a,b)=> meta[b].peak_height - meta[a].peak_height);
    } else if(mode==='peak_time'){
      idx.sort((a,b)=> meta[b].peak_time - meta[a].peak_time);
    } else if(mode==='closest_to_ref'){
      idx.sort((a,b)=>{
        var da = Math.abs(meta[a].peak_time - t_ref);
        var db = Math.abs(meta[b].peak_time - t_ref);
        if(da !== db) return da - db;
        return meta[b].peak_height - meta[a].peak_height;
      });
    }
    return idx;
  }

  // Only does ranking + visibility + table. NO relayout here.
  function updateRankingAndVisibility(tref){
    var mode = document.getElementById('pi_sort_mode').value;
    var topn = parseInt(document.getElementById('pi_top_n').value);

    var ord = order_indices(mode, tref);

    var vis = new Array(meta.length).fill(false);
    for(var k=0; k<Math.min(topn, ord.length); k++){
      vis[ord[k]] = true;
    }
    Plotly.restyle(gd, 'visible', vis);

    var rows = ord.slice(0, Math.min(topn, ord.length)).map(function(i, rnk){
      var d = meta[i];
      var peaked = d.peak_time < tref;
      return '<tr>'
        + '<td>'+(rnk+1)+'</td>'
        + '<td>'+d.locus+'</td>'
        + '<td>'+d.n_sites+'</td>'
        + '<td>'+d.peak_height.toPrecision(6)+'</td>'
        + '<td>'+d.peak_time.toPrecision(6)+'</td>'
        + '<td>'+(peaked ? 'TRUE' : 'FALSE')+'</td>'
        + '</tr>';
    }).join('');

    document.getElementById('pi_rank_body').innerHTML = rows;
    document.getElementById('pi_top_n_label').innerText = topn.toString();
    document.getElementById('pi_t_ref_label').innerText = tref.toPrecision(6);
  }

  // Slider-driven: update the dashed line, then update ranking/visibility.
  function apply_controls(){
    var tref = parseFloat(document.getElementById('pi_t_ref').value);

    // Update vertical reference line (this will emit plotly_relayout)
    Plotly.relayout(gd, {'shapes[0].x0': tref, 'shapes[0].x1': tref});

    // Update ranking/visibility/table WITHOUT further relayout
    updateRankingAndVisibility(tref);
  }

  // UI events
  document.getElementById('pi_sort_mode').addEventListener('change', function(){
    var tref = parseFloat(document.getElementById('pi_t_ref').value);
    updateRankingAndVisibility(tref);
  });
  document.getElementById('pi_top_n').addEventListener('input', function(){
    var tref = parseFloat(document.getElementById('pi_t_ref').value);
    updateRankingAndVisibility(tref);
  });
  document.getElementById('pi_t_ref').addEventListener('input', apply_controls);

  // Export buttons (Route A: no Python dependency)
  // Guard against missing buttons so the HTML never fails to load.
  var btnSvg = document.getElementById('pi_download_svg');
  if(btnSvg){
    btnSvg.addEventListener('click', function(){
      Plotly.downloadImage(gd, {format:'svg', filename:'PI_multi_locus', scale:2});
    });
  }

  var btnPdf = document.getElementById('pi_print_pdf');
  if(btnPdf){
    btnPdf.addEventListener('click', function(){
      window.print(); // user selects 'Save as PDF'
    });
  }

  // Line-drag driven: relayout event tells us new x0, we sync slider + update ranking.
  // NOTE: This requires plotly config editable=TRUE on the R side.
  gd.on('plotly_relayout', function(ev){
    if(ev['shapes[0].x0'] !== undefined){
      var tref = parseFloat(ev['shapes[0].x0']);

      // Sync slider value (setting .value does NOT fire input event automatically)
      document.getElementById('pi_t_ref').value = tref;

      // Update ranking/visibility/table WITHOUT relayout
      updateRankingAndVisibility(tref);

      // Keep label consistent
      document.getElementById('pi_t_ref_label').innerText = tref.toPrecision(6);
    }
  });

  // Initial render
  apply_controls();
}
", jsonlite::toJSON(meta, dataframe = "rows", auto_unbox = TRUE))
  
  
  ## UI controls (adds SVG export + browser print-to-PDF)
  controls <- htmltools::tags$div(
    class = "pi-controls",
    style = paste(
      "font-family: Arial, sans-serif;",
      "margin: 10px 0;",
      "position: sticky; top: 0; z-index: 9999;",   # keep controls visible while scrolling
      "background: white; padding: 8px 6px;",
      "border-bottom: 1px solid #ddd;",
      sep = " "
    ),
    htmltools::tags$div(
      style = "display:flex; gap:18px; align-items:center; flex-wrap: wrap;",
      
      ## Sort mode
      htmltools::tags$label(
        style = "display:flex; gap:8px; align-items:center;",
        "Sort loci by:",
        htmltools::tags$select(
          id = "pi_sort_mode",
          htmltools::tags$option(value = "peak_height", "Peak height (desc)"),
          htmltools::tags$option(value = "peak_time", "Peak time (latest first)"),
          htmltools::tags$option(value = "closest_to_ref", "Peak closest to selected time")
        )
      ),
      
      ## Top-N slider
      htmltools::tags$label(
        style = "display:flex; gap:8px; align-items:center;",
        "Show top N:",
        htmltools::tags$input(
          id = "pi_top_n", type = "range",
          min = 1, max = length(loci), value = show_n, step = 1,
          style = "width: 220px;"
        ),
        htmltools::tags$span(id = "pi_top_n_label", show_n)
      ),
      
      ## Time-of-interest slider (moves dashed line + re-ranks)
      htmltools::tags$label(
        style = "display:flex; gap:8px; align-items:center;",
        "Time of interest:",
        htmltools::tags$input(
          id = "pi_t_ref", type = "range",
          min = min(times), max = max(times), value = t_ref0,
          step = (max(times) - min(times)) / 200,
          style = "width: 260px;"
        ),
        htmltools::tags$span(id = "pi_t_ref_label", signif(t_ref0, 6))
      ),
      
      ## Export buttons (Route A: no Python/reticulate dependency)
      htmltools::tags$div(
        style = "display:flex; gap:10px; align-items:center;",
        htmltools::tags$button(
          id = "pi_download_svg",
          type = "button",
          style = "padding:6px 10px; cursor:pointer;",
          "Download SVG (publication)"
        ),
        htmltools::tags$button(
          id = "pi_print_pdf",
          type = "button",
          style = "padding:6px 10px; cursor:pointer;",
          "Save as PDF (browser print)"
        )
      )
    )
  )
  
  
  ranking_block <- htmltools::tags$div(
    class = "pi-section",
    htmltools::tags$b("Current ranking (top N):"),
    htmltools::tags$table(
      style="border-collapse: collapse; width: 100%; margin-top: 8px;",
      htmltools::tags$thead(
        htmltools::tags$tr(
          lapply(c("Rank","Locus","n_sites","peak_height","peak_time","peaked_by_ref"),
                 function(h) htmltools::tags$th(h, style="border:1px solid #ddd; padding:6px; text-align:left;"))
        )
      ),
      htmltools::tags$tbody(id="pi_rank_body")
    )
  )
  
  ## Build a single htmlwidget (required by saveWidget)
  p_widget <- htmlwidgets::onRender(p, js)
  
  ## Inject CSS + controls + ranking block into the widget's HTML
  p_widget <- htmlwidgets::prependContent(p_widget, style_tag, controls, ranking_block)
  
  ## Save as a standalone HTML file
  if (!is.null(file)) {
    htmlwidgets::saveWidget(p_widget, file = file, selfcontained = selfcontained)
  }
  
  return(p_widget)
}



## ---------------------------------------------------------------------------
## Interactive tree + per-edge multi-locus signal-minus-noise lollipop (HTML)
## Interpretation 2A (Route A export):
## - User clicks an edge midpoint on the phylogeny (left panel)
## - Right panel updates a lollipop plot across loci:
##     baseline at 0, stem up/down to (S-N), dot at (S-N)
##     where S-N = P_correct - P_wrong (same math as legacy Approximator)
## - Export: client-side SVG download (user can print SVG to PDF)
##
## NOTES:
## - Algorithm unchanged: uses psnr() and pnl() from the original code.
## - Exact speedup: folds rates into unique values + multiplicities (no approximation).
## - UI: scroll-zoom disabled to avoid trackpad/wheel mis-zoom; focused on click + controls.
## ---------------------------------------------------------------------------

## ---- small helper (only once; safe if already defined) ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

## ---- fold rates (exact) ----
.fold_rates <- function(rate_vec) {
  rate_vec <- as.numeric(rate_vec)
  tab <- table(rate_vec)
  list(u = as.numeric(names(tab)), c = as.numeric(tab))
}

## ---- folded Approximator (exact; same algebra as legacy Approximator) ----
Approximator_folded <- function(t, t0, u_rates, u_counts, s) {
  if (!is.finite(t) || !is.finite(t0) || t < 0 || t0 < 0) {
    return(c(P_correct = NA_real_, P_poly = NA_real_, P_wrong = NA_real_))
  }
  
  ## IMPORTANT: psnr() and pnl() must already exist in PhyInformR.R
  npsnr <- psnr(u_rates, t, t0, s)
  npnl  <- pnl(u_rates,  t, t0, s)
  
  ## Folded sums (exact)
  eYsum   <- sum(u_counts * npsnr)
  eX1sum  <- sum(u_counts * npnl)
  eY2sum  <- sum(u_counts * (npsnr^2))
  eX12sum <- sum(u_counts * (npnl^2))
  eX1Ysum <- sum(u_counts * (npsnr * npnl))
  eSQsum  <- sum(u_counts * (npsnr * sqrt(pmax(npnl, 0))))  ## numeric guard only
  
  ## Same algebra as legacy Approximator
  Expectationx <- eX1sum + sqrt(eX1sum / pi)
  Expectation  <- eYsum - Expectationx
  
  variancey <- eYsum - eY2sum
  variancex <- ((pi - 1) / pi) * eX1sum - eX12sum
  variance  <- variancey + variancex - 2 * eX1Ysum - (2 / sqrt(pi)) * eSQsum
  
  ## Guard: variance can be tiny negative from floating error
  variance <- pmax(variance, 0)
  sdv <- sqrt(variance)
  
  pr_wrong <- stats::pnorm(-0.5, mean = Expectation, sd = sdv)
  pr_poly  <- stats::pnorm( 0.5, mean = Expectation, sd = sdv) -
    stats::pnorm(-0.5, mean = Expectation, sd = sdv)
  pr_corr  <- 1 - stats::pnorm(0.5, mean = Expectation, sd = sdv)
  
  c(P_correct = as.numeric(pr_corr),
    P_poly    = as.numeric(pr_poly),
    P_wrong   = as.numeric(pr_wrong))
}

## ---------------------------------------------------------------------------
## Main entry: interactive tree + edge click -> per-locus S-N lollipop plot (HTML)
## ---------------------------------------------------------------------------
tree_signal_noise_multi_html <- function(rates_list, tree,
                                         s = 4,
                                         file = NULL,
                                         selfcontained = TRUE,
                                         default_top_n = 25) {
  
  if (!requireNamespace("ape", quietly = TRUE)) stop("Please install.packages('ape')")
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Please install.packages('plotly')")
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) stop("Please install.packages('htmlwidgets')")
  if (!requireNamespace("htmltools", quietly = TRUE)) stop("Please install.packages('htmltools')")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Please install.packages('jsonlite')")
  
  if (!is.list(rates_list) || length(rates_list) < 2) {
    stop("rates_list must be a named list with >=2 loci (genes/partitions).")
  }
  if (is.null(names(rates_list)) || any(names(rates_list) == "")) {
    names(rates_list) <- paste0("Locus_", seq_along(rates_list))
  }
  
  ## ---- Node/edge metadata (internal edges only, per your NodeWalker) ----
  if (!exists("NodeWalker")) stop("NodeWalker(tree) not found in the environment. Ensure PhyInformR.R defines NodeWalker().")
  nodes <- NodeWalker(tree)  ## must include: parent_node, daughter_node, p_node_time, d_node_time
  
  req_cols <- c("parent_node", "daughter_node", "p_node_time", "d_node_time")
  miss_cols <- setdiff(req_cols, colnames(nodes))
  if (length(miss_cols) > 0) {
    stop("NodeWalker(tree) is missing required columns: ", paste(miss_cols, collapse = ", "))
  }
  
  n_edges <- nrow(nodes)
  loci    <- names(rates_list)
  n_loci  <- length(loci)
  
  ## ---- Pre-fold all loci once (exact, faster; algorithm unchanged) ----
  folded <- lapply(rates_list, .fold_rates)
  
  ## ---- Compute P_correct, P_wrong, and S-N matrices: edges x loci ----
  ## Edge definition follows legacy PlotTreeSI: t = daughter_time, t0 = parent_time - daughter_time.
  pc_mat <- matrix(NA_real_, nrow = n_edges, ncol = n_loci, dimnames = list(NULL, loci))
  pw_mat <- matrix(NA_real_, nrow = n_edges, ncol = n_loci, dimnames = list(NULL, loci))
  sn_mat <- matrix(NA_real_, nrow = n_edges, ncol = n_loci, dimnames = list(NULL, loci))
  
  for (j in seq_len(n_loci)) {
    fu <- folded[[j]]$u
    fc <- folded[[j]]$c
    for (i in seq_len(n_edges)) {
      t  <- as.numeric(nodes[i, "d_node_time"])
      t0 <- as.numeric(nodes[i, "p_node_time"]) - t
      probs <- Approximator_folded(t, t0, fu, fc, s)
      pc <- as.numeric(probs["P_correct"])
      pw <- as.numeric(probs["P_wrong"])
      pc_mat[i, j] <- pc
      pw_mat[i, j] <- pw
      sn_mat[i, j] <- pc - pw
    }
  }
  
  ## ---- Tree coordinates (no device output) ----
  ape::plot.phylo(tree, plot = FALSE, show.tip.label = FALSE, direction = "leftwards")
  PlotEnv <- get(".PlotPhyloEnv", envir = asNamespace("ape"))
  lp <- get("last_plot.phylo", envir = PlotEnv)
  xx <- lp$xx
  yy <- lp$yy
  
  ## ---- Build edge segments for drawing ----
  ed <- tree$edge
  seg_x <- c()
  seg_y <- c()
  for (k in seq_len(nrow(ed))) {
    p <- ed[k, 1]
    d <- ed[k, 2]
    seg_x <- c(seg_x, xx[p], xx[d], NA)
    seg_y <- c(seg_y, yy[p], yy[d], NA)
  }
  
  ## ---- Build clickable edge midpoints (mapped by NodeWalker parent/daughter) ----
  mid_x <- numeric(n_edges)
  mid_y <- numeric(n_edges)
  for (i in seq_len(n_edges)) {
    p <- as.integer(nodes[i, "parent_node"])
    d <- as.integer(nodes[i, "daughter_node"])
    mid_x[i] <- (xx[p] + xx[d]) / 2
    mid_y[i] <- (yy[p] + yy[d]) / 2
  }
  
  edge_meta <- data.frame(
    edge_index    = seq_len(n_edges),
    parent_node   = as.integer(nodes[, "parent_node"]),
    daughter_node = as.integer(nodes[, "daughter_node"]),
    p_time        = as.numeric(nodes[, "p_node_time"]),
    d_time        = as.numeric(nodes[, "d_node_time"]),
    x             = mid_x,
    y             = mid_y
  )
  
  ## ---- Initial selection (edge 1) ----
  edge0 <- 1L
  sn0 <- as.numeric(sn_mat[edge0, ])
  pc0 <- as.numeric(pc_mat[edge0, ])
  pw0 <- as.numeric(pw_mat[edge0, ])
  
  ## ---- Helper: build top-N payload ----
  make_top_payload <- function(sn, pc, pw, loci, top_n = default_top_n, sort_mode = "abs_desc") {
    df <- data.frame(locus = loci, S = pc, N = pw, SN = sn, stringsAsFactors = FALSE)
    if (sort_mode == "sn_desc") {
      df <- df[order(df$SN, decreasing = TRUE), ]
    } else if (sort_mode == "sn_asc") {
      df <- df[order(df$SN, decreasing = FALSE), ]
    } else { ## abs_desc
      df <- df[order(abs(df$SN), decreasing = TRUE), ]
    }
    top_n <- min(top_n, nrow(df))
    df <- df[seq_len(top_n), ]
    df$rank <- seq_len(nrow(df))
    df
  }
  
  default_top_n <- min(as.integer(default_top_n), n_loci)
  top0 <- make_top_payload(sn0, pc0, pw0, loci, top_n = default_top_n, sort_mode = "abs_desc")
  
  ## Lollipop stems as x,y with NA breaks
  stem_x0 <- as.vector(rbind(top0$rank, top0$rank, rep(NA, nrow(top0))))
  stem_y0 <- as.vector(rbind(rep(0, nrow(top0)), top0$SN, rep(NA, nrow(top0))))
  
  ## ---- Build plotly: two-panel subplot (tree | lollipop) ----
  ## Panel A: tree
  p_tree <- plotly::plot_ly() %>%
    plotly::add_trace(
      x = seg_x, y = seg_y,
      type = "scatter", mode = "lines",
      line = list(width = 1),
      hoverinfo = "skip",
      name = "tree"
    ) %>%
    plotly::add_trace(
      data = edge_meta,
      x = ~x, y = ~y,
      type = "scatter", mode = "markers",
      marker = list(size = 7, opacity = 0.85),
      customdata = ~edge_index,
      text = ~paste0(
        "edge=", edge_index,
        "<br>parent_time=", signif(p_time, 6),
        "<br>daughter_time=", signif(d_time, 6)
      ),
      hoverinfo = "text",
      name = "click edge"
    ) %>%
    plotly::layout(
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      margin = list(l = 10, r = 10, t = 40, b = 10),
      title = list(text = "Phylogeny (click an edge midpoint)")
    )
  
  ## Panel B: lollipop (baseline at 0, stems + dots)
  p_lol <- plotly::plot_ly() %>%
    plotly::add_trace(
      x = stem_x0, y = stem_y0,
      type = "scatter", mode = "lines",
      hoverinfo = "skip",
      name = "stem"
    ) %>%
    plotly::add_trace(
      data = top0,
      x = ~rank, y = ~SN,
      type = "scatter", mode = "markers",
      text = ~paste0(
        "locus=", locus,
        "<br>P_correct(S)=", signif(S, 6),
        "<br>P_wrong(N)=", signif(N, 6),
        "<br>S-N=", signif(SN, 6)
      ),
      hoverinfo = "text",
      name = "S-N"
    ) %>%
    plotly::layout(
      xaxis = list(title = "Locus (ranked)"),
      yaxis = list(title = "Signal − Noise (S − N)"),
      shapes = list(list(
        type = "line",
        x0 = 0.5, x1 = nrow(top0) + 0.5,
        y0 = 0,   y1 = 0,
        xref = "x", yref = "y",
        line = list(width = 2)
      )),
      margin = list(l = 70, r = 10, t = 40, b = 55)
    )
  
  ## Subplot
  p <- plotly::subplot(p_tree, p_lol,
                       widths = c(0.45, 0.55),
                       shareY = FALSE, titleX = TRUE, titleY = TRUE)
  
  ## ---- Controls (sorting + top N + export) ----
  controls <- htmltools::tags$div(
    class = "sn-controls",
    style = "font-family: Arial, sans-serif;",
    htmltools::tags$div(
      style = "display:flex; gap:20px; align-items:center; flex-wrap: wrap;",
      htmltools::tags$label("Sort loci by: ",
                            htmltools::tags$select(
                              id = "sn_sort_mode",
                              htmltools::tags$option(value = "abs_desc", "abs(S−N) (desc)"),
                              htmltools::tags$option(value = "sn_desc", "S−N (desc)"),
                              htmltools::tags$option(value = "sn_asc",  "S−N (asc)")
                            )
      ),
      htmltools::tags$label("Show top N: ",
                            htmltools::tags$input(
                              id = "sn_top_n", type = "range",
                              min = 1, max = n_loci, value = default_top_n, step = 1
                            ),
                            htmltools::tags$span(id = "sn_top_n_label", default_top_n)
      ),
      htmltools::tags$button(
        id = "sn_export_svg",
        type = "button",
        "Export SVG (then print to PDF)"
      ),
      htmltools::tags$span(id = "sn_edge_label", style = "margin-left:10px;")
    )
  )
  
  ## Page-level CSS (scroll + sticky header)
  style_tag <- htmltools::tags$style(htmltools::HTML("
    html, body { height: auto; overflow-y: auto; margin: 0; padding: 0; }
    .sn-controls {
      position: sticky;
      top: 0;
      background: #ffffff;
      z-index: 9999;
      border-bottom: 1px solid #e5e5e5;
      padding: 10px 12px;
    }
    .sn-controls button { padding: 6px 10px; }
    .sn-controls label { margin-right: 12px; }
    .sn-controls input[type=range] { vertical-align: middle; }
  "))
  
  ## ---- JS: click edge -> update lollipop (stems + dots) ----
  ## Trace indices inside subplot:
  ## 0: tree lines, 1: tree click markers, 2: lollipop stems, 3: lollipop dots
  ##
  ## IMPORTANT:
  ## - Avoid relayout loops: only update baseline length + annotation text.
  ## - Scroll zoom disabled via config.
  js <- sprintf("
  function(el, x){
    var gd = document.getElementById(el.id);

    var loci  = %s;
    var pcMat = %s;
    var pwMat = %s;
    var snMat = %s;

    var currentEdge = 1;

    function makeTop(edgeIdx){
      var mode = document.getElementById('sn_sort_mode').value;
      var topn = parseInt(document.getElementById('sn_top_n').value);

      var sn = snMat[edgeIdx-1];
      var pc = pcMat[edgeIdx-1];
      var pw = pwMat[edgeIdx-1];

      var arr = loci.map(function(name, i){
        return {locus:name, S:pc[i], N:pw[i], SN:sn[i]};
      });

      if(mode === 'sn_desc'){
        arr.sort(function(a,b){ return b.SN - a.SN; });
      } else if(mode === 'sn_asc'){
        arr.sort(function(a,b){ return a.SN - b.SN; });
      } else {
        arr.sort(function(a,b){ return Math.abs(b.SN) - Math.abs(a.SN); });
      }

      topn = Math.min(topn, arr.length);
      arr = arr.slice(0, topn);
      arr.forEach(function(d,i){ d.rank = i+1; });

      return arr;
    }

    function ensureAnno(){
      var hasAnno = (gd.layout.annotations && gd.layout.annotations.length > 0);
      if(!hasAnno){
        Plotly.relayout(gd, {
          annotations: [{
            xref: 'paper', yref: 'paper',
            x: 0.73, y: 1.05,
            showarrow: false,
            text: 'Per-locus S−N (edge 1)',
            font: {size: 14}
          }]
        });
      }
    }

    function updateLollipop(edgeIdx){
      currentEdge = edgeIdx;

      var top = makeTop(edgeIdx);

      var stemX = [];
      var stemY = [];
      for(var i=0; i<top.length; i++){
        stemX.push(top[i].rank); stemY.push(0);
        stemX.push(top[i].rank); stemY.push(top[i].SN);
        stemX.push(null);        stemY.push(null);
      }

      var dotX = top.map(function(d){ return d.rank; });
      var dotY = top.map(function(d){ return d.SN; });
      var dotText = top.map(function(d){
        return 'locus=' + d.locus
          + '<br>P_correct(S)=' + Number(d.S).toPrecision(6)
          + '<br>P_wrong(N)=' + Number(d.N).toPrecision(6)
          + '<br>S-N=' + Number(d.SN).toPrecision(6);
      });

      // Update traces 2 and 3
      Plotly.restyle(gd, {x:[stemX], y:[stemY]}, [2]);
      Plotly.restyle(gd, {x:[dotX],  y:[dotY],  text:[dotText]}, [3]);

      // Update baseline length and title annotation (no zoom/pan changes)
      Plotly.relayout(gd, {
        'shapes[0].x0': 0.5,
        'shapes[0].x1': (top.length + 0.5),
        'annotations[0].text': ('Per-locus S−N (edge ' + edgeIdx + ')')
      });

      document.getElementById('sn_top_n_label').innerText = top.length.toString();
      document.getElementById('sn_edge_label').innerText  = 'Selected edge: ' + edgeIdx;
    }

    // Handle clicks: only respond to edge midpoint markers (trace 1)
    gd.on('plotly_click', function(ev){
      if(!ev || !ev.points || ev.points.length === 0) return;
      var pt = ev.points[0];
      if(pt.curveNumber !== 1) return;
      var edgeIdx = pt.customdata;
      if(edgeIdx === undefined || edgeIdx === null) return;
      updateLollipop(parseInt(edgeIdx));
    });

    // UI events
    document.getElementById('sn_sort_mode').addEventListener('change', function(){
      updateLollipop(currentEdge);
    });
    document.getElementById('sn_top_n').addEventListener('input', function(){
      updateLollipop(currentEdge);
    });

    // Export SVG (Route A)
    document.getElementById('sn_export_svg').addEventListener('click', function(){
      Plotly.downloadImage(gd, {format:'svg', filename:'tree_signal_noise'});
    });

    // Initialize
    ensureAnno();
    document.getElementById('sn_edge_label').innerText = 'Selected edge: 1';
    updateLollipop(1);
  }
  ",
                jsonlite::toJSON(loci, auto_unbox = TRUE),
                jsonlite::toJSON(unname(pc_mat), dataframe = "values", auto_unbox = TRUE),
                jsonlite::toJSON(unname(pw_mat), dataframe = "values", auto_unbox = TRUE),
                jsonlite::toJSON(unname(sn_mat), dataframe = "values", auto_unbox = TRUE)
  )
  
  p_widget <- htmlwidgets::onRender(p, js)
  p_widget <- htmlwidgets::prependContent(p_widget, style_tag, controls)
  
  ## ---- Plotly config: avoid accidental zoom; keep focused interactions ----
  ## NOTE: do NOT set editable=TRUE here (causes drag handles / shape editing).
  p_widget$x$config <- modifyList(
    p_widget$x$config %||% list(),
    list(
      scrollZoom = FALSE,          ## critical: disable mouse-wheel/trackpad zoom
      doubleClick = "reset",
      displaylogo = FALSE,
      editable = FALSE,
      modeBarButtonsToRemove = c(
        "zoom2d","select2d","lasso2d",
        "autoScale2d","resetScale2d",
        "zoomIn2d","zoomOut2d","pan2d"
      )
    )
  )
  
  if (!is.null(file)) {
    htmlwidgets::saveWidget(p_widget, file = file, selfcontained = selfcontained)
  }
  return(p_widget)
}
