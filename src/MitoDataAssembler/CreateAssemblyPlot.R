#bounds of circle
bounds=c(-100,100)
#outerRadius where text annotation starts
outerRadius=90
innerRadius=90
genomeSize=16569

createLayout<-function(titName,genomeSize)
{
  #plot(bounds,bounds,type="n",yaxt="n",xaxt="n",ann=FALSE)
  par(mar=c(0,0,2,0))
  plot.new()
  plot.window(xlim=bounds,ylim=bounds)
  box() 
  radians=seq(0,2*pi,0.01)
  x=innerRadius*cos(radians)
  y=innerRadius*sin(radians)
  #lines(x,y)
  polygon(x,y)  
  #add ticks and annotations every
  interval=2500
  pos=seq(0,genomeSize,interval)
  posr=pi/2-2*pi*pos/genomeSize  
  x1=innerRadius*cos(posr)
  y1=innerRadius*sin(posr)
  textRadius=innerRadius+3
  x2=textRadius*cos(posr)
  y2=textRadius*sin(posr)
  x=cbind(x1,x2)
  y=cbind(y1,y2)
  #now add lines 
  for(i in 1:length(x1))
  {
    lines(x[i,],y[i,])
   
  }
  textRadius=innerRadius+9
  x2=textRadius*cos(posr)
  y2=textRadius*sin(posr)
  #now add lines 
  for(i in 1:length(x2))
  {
    text(x2[i],y2[i],pos[i])
  }    
  title(titName)
}

createLayout("Genome Assembly")


t