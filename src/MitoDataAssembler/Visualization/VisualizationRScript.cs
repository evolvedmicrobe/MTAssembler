﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MitoDataAssembler.Visualization
{
    public static class VisualizationRScript
    {
        /// <summary>
        /// A cut and paste version of the RSCript with '"' changed to '""'.  Always directly edit the R Script.
        /// </summary>
        public static string FILE_AS_STRING = @"
##Script to create a genome assembly from the whole thing
library(grDevices)
library(grid)
####GENERAL PARAMETERS
filename=""test.pdf""
assemblyLength=16569
#This value is hard coded in c#
outerRadius=.45
init.angle=pi/2
radian.delta=10;
radius.delta=1/1000;
tickInterval=2500
minAngleForStraightArrow=pi/12
arrowLoopDif=0.03
arrowDeltaIntervals=100
###GENERAL FUNCTIONS

####INITIAL LAYOUT
getTickCoords<-function(x,mags) {  
  t(apply(rbind(mags,rep(x,3)),2,polar2xy))
}

drawTick<-function(x,mags) {
  xys=getTickCoords(x,mags)
  xys=xys
  grid.lines(x=xys[1:2,1],y=xys[1:2,2])
  val=paste(format(x,big.mark="",""),""bp"",sep="" "")
  grid.text(val,x=xys[3,1],y=xys[3,2])
}

#takes a magnitude/angle
polar2xy <- function(vals) {
  magnitude=vals[1]
  dAngle=vals[2]
  #get angle
  ang <- init.angle - (dAngle/assemblyLength)* 2* pi
  #convert to x,y
  #list(x = magnitude* cos(ang), y = magnitude * sin(ang))
  xy=.5+c(magnitude*cos(ang), magnitude*sin(ang))
  return(xy)
}
createPlot<-function(filename)
{
  #windows()  
  
  pdf(filename,paper=""letter"")
  grid.newpage()
  #create a two panel layout, forcing the middle viewport to be square
  vlayout<-grid.layout(nrow=2,ncol=1,heights=unit(c(2,1),c(""lines"",""null"")),respect=c(0,1))
  vp1<-viewport(layout=vlayout)
  pushViewport(vp1)
  #add title bar
  vp2<-viewport(layout.pos.row=1,name=""Title"")
  pushViewport(vp2)
  grid.text(""Mitochondrial Genome Assembly"")
  popViewport()
  #Get Settings
  gps<-get.gpar()
  #change font size and fill
  gps$fontsize<-12
  gps$fill<-""transparent""
  #make new viewport
  vp3<-viewport(layout.pos.row=2,name=""Assembly"",gp=gps)
  pushViewport(vp3)
  grid.circle(r=outerRadius)
  
  #### ADD TICK MARKS
  ticks=seq(0,assemblyLength,tickInterval) 
  endTick=outerRadius*1.02;
  tickText=outerRadius*1.1
  mags=c(outerRadius,endTick,tickText)
  tickDraw<-function(x){ drawTick(x,mags)}
  sapply(ticks,tickDraw)
  
  
}

###PLOT ASSEMBLED REGIONS
createArc<-function(magnitude,sdAngle,edAngle) {
  if(edAngle<sdAngle){
    angs=c(seq(sdAngle,assemblyLength,radian.delta),seq(0,edAngle,radian.delta))
  }
  else {angs=seq(sdAngle,edAngle,radian.delta)}
  toChange=cbind(rep(magnitude,length(angs)),angs)
  coords=t(apply(toChange,1,polar2xy))
  #grid.lines(coords[,1],coords[,2])
}

createRay<-function(ang,s,e) {
  if(s>e) {ray=seq(s,e,-radius.delta)} else 
  {ray=seq(s,e,radius.delta)}
  toChange=cbind(ray,rep(ang,length(ray)))
  coords=t(apply(toChange,1,polar2xy))
}




drawSegment<-function(start,end,low,height) {
  #print(c(start,end,low,height))
  fill=get.gpar()
  fill$fill=""blue""
  fill$lwd=2
  ray1=createRay(start,low,low+height)
  arc1=createArc(low+height,start,end)
  ray2=createRay(end,low+height,low)
  arc2=createArc(low,start,end)
  arc2=arc2[seq(nrow(arc2),1,-1),]
  path=rbind(ray1,arc1,ray2,arc2)
  grid.polygon(path[,1],path[,2],gp=fill)
}

reduceAngle<-function(ang) {
  ang=2*pi*(ang/assemblyLength)
  ang=ang-2*pi*floor(ang/(2*pi))
  return(ang)
}

#start of function to draw as loops, opted to give it up though
drawLoopedArrow<-function(sMag,sAng,eMag,eAng,arrowEnd) {
  #This is a looped arrow that will go out, then back in
  #angles from and two
  
  #get coordinates for a point above the current point
  r=.5*arrowLoopDif
  np=sMag+r
  na=sAng
  xy=polar2xy(c(np,na))
  degO=360*sAng/assemblyLength
  #Convert
  ndeg=360-degO
  #now make a viewpoint at this location
  vpc=viewport(x=xy[1],y=xy[2],width=arrowLoopDif,height=arrowLoopDif,just=c(""center"",""center""),angle=ndeg)
  pushViewport(vpc)
  #now we have the viewport
  arc1=createArc(1,assemblyLength/2,assemblyLength)
  grid.lines(arc1[,1],arc1[,2])
  #now where is the end value?
  
  xy2=polar2xy(c(eMag,eAng))
  deltxy=xy2-xy
  #get hypotenus and square
  mag=sqrt(sum(deltxy^2))*arrowLoopDif
  ang=atan(deltxy[2]/deltxy[1])
  ang=2*pi-ang+pi/4
  #recalibrate angle
  drawStraightArrow(1,0,mag,ang,TRUE)
  popViewport()
}


drawStraightArrow<-function(sMag,sAng,eMag,eAng,arrowAtEnd) {
  d1=(eAng-sAng)/arrowDeltaIntervals
  angs=seq(sAng,eAng,d1)
  d1=(eMag-sMag)/arrowDeltaIntervals
  mags=seq(sMag,eMag,d1)
  #angs=c(angs,eAng) #this screwed up arrow head for some reason
  #mags=c(mags,eMag)
  toChange=cbind(mags,angs)
  if(!arrowAtEnd) {toChange=toChange[seq(nrow(toChange),1,-1),]}
  coords=t(apply(toChange,1,polar2xy))
  arr=arrow()
  #grid.arrows(coords[,1],coords[,2])
  grid.lines(coords[,1],coords[,2])
  n=length(angs)
  last2=(n-1):n
  gp=get.gpar()
  gp$lwd=4
  gp$col=""red""
  grid.lines(coords[last2,1],coords[last2,2],arrow=arr,gp=gp)
}

drawArrow<-function(sMag,sAng,eMag,eAng) {
  #Two options here, if close we loop the arrow, if far, 
  #direct path, so calculate distance first   
  #get reduced angles on 0-2pi scale
  sAng=reduceAngle(sAng)
  eAng=reduceAngle(eAng)
  arrowAtEnd=TRUE
  #now orient so ang2 is largest
  if(sAng>eAng){
    tmp=c(sMag,sAng)
    sMag=eMag;sAng=eAng;eAng=tmp[2];eMag=tmp[1];
    arrowAtEnd=FALSE
  }
  #now do we go direct or loop around?
  opt1=eAng-sAng #direct path
  opt2=sAng+(2*pi)-eAng #loop path
  if(opt2<opt1) {
    #flip again
    sAng=sAng+2*pi
    tmp=c(sMag,sAng)
    sMag=eMag;sAng=eAng;eAng=tmp[2];eMag=tmp[1];
    arrowAtEnd=!arrowAtEnd
    }
  
  #now decide which way we draw, 
  
  if(abs(eAng-sAng)<minAngleForStraightArrow) {
    arrowFunc=drawLoopedArrow
  } else {
    arrowFunc=drawStraightArrow
  }
  convF=assemblyLength/(2*pi)
  arrowFunc(sMag,sAng*convF,eMag,eAng*convF,arrowAtEnd)
}
#grid.newpage()
#drawArrow(.25,15000,.1,2500)
#drawStraightArrow(.2,15000,.4,15001,TRUE)
#drawSegment(6000,15000,x,.4-x)

demoDrawing<-function()
{
  createPlot(""test.pdf"")
  nsegs=rbinom(1,10,.5)
  inRad=.1
  out=outerRadius
  height=(out-inRad)/nsegs
  positions=inRad+(0:(nsegs-1))*height
  size=16569
  angs=c()
  for(start in positions)
  {
    ang=runif(2)*16569
    ang=c(min(ang),max(ang))
    drawSegment(ang[1],ang[2],start,height-height*.1)
    angs=rbind(angs,ang)
  }
  allPoints=1:nsegs
  #now random arrows
  for(r in 1:4)
  {
    startEnd=sample(allPoints,2)
    f=startEnd[1]
    e=startEnd[2]
    start=angs[f,]
    end=angs[e,]
    fh=inRad+(f-1)*height+(height/2)
    eh=inRad+(e-1)*height+(height/2)
    #sMag,sAng,eMag,eAng
    drawArrow(fh,angs[f,1],eh,angs[e,2])    
  }
  #dev.off()
}    
#demoDrawing()





#dev.off()
        ";
    }
}
