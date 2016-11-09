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
layout(mat=mat,heights=c(250,300))
par(mar=c(0,0,0,0), oma=c(5,5,1,1))
#par(bg = "white")   
#split.screen(c(2,1))
#screen(1)
#plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
##coord are left,right,bottom,top from 0 to 1

par(plt=c(0,0.9,0.2,0.99))
plot(tree,show.tip.label=FALSE,direction="l")
####Lower corner, note that the pi is offset to mirror the trees end

par(plt=c(0.027,0.9,0,0.99))

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
layout(mat=mat,heights=c(250,300))
par(mar=c(0,0,0,0), oma=c(5,5,1,1))
#par(bg = "white")   
#split.screen(c(2,1))
#screen(1)
#plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
##coord are left,right,bottom,top from 0 to 1

par(plt=c(0,0.9,0.2,0.99))
plot(tree,show.tip.label=FALSE,direction="l")
####Lower corner, note that the pi is offset to mirror the trees end

par(plt=c(0.027,0.9,0,0.99))

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
layout(mat=mat,heights=c(250,300))
par(mar=c(0,0,0,0), oma=c(5,5,1,1))
#par(bg = "white")   
#split.screen(c(2,1))
#screen(1)
#plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
##coord are left,right,bottom,top from 0 to 1

par(plt=c(0,0.9,0.2,0.99))
plot(tree,show.tip.label=FALSE,direction="l")
####Lower corner, note that the pi is offset to mirror the trees end

par(plt=c(0.027,0.9,0,0.99))

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
layout(mat=mat,heights=c(250,300))
par(mar=c(0,0,0,0), oma=c(5,5,1,1))
#par(bg = "white")   
#split.screen(c(2,1))
#screen(1)
#plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
##coord are left,right,bottom,top from 0 to 1

par(plt=c(0,0.9,0.2,0.99))
plot(tree,show.tip.label=FALSE,direction="l")
####Lower corner, note that the pi is offset to mirror the trees end

par(plt=c(0.027,0.9,0,0.99))

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
layout(mat=mat,heights=c(250,300))
par(mar=c(0,0,0,0), oma=c(5,5,1,1))
par(plt=c(0,0.9,0.2,0.99))
plot(tree,show.tip.label=FALSE,direction="l")
par(plt=c(0.027,0.9,0,0.99))

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
	par(plt=c(0.027,0.9,0,0.99))

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
	rbind(princtree$value, prpolytomy$value, prcortree$value)->output
	return(output)
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
  print(qp2)
  print(stored_ints)
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
