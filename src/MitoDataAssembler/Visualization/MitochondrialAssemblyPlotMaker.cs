using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly.Padena;
using Bio.Algorithms.Assembly.Graph;
using MitoDataAssembler.Extensions;
using MitoDataAssembler.R;
using System.Diagnostics;
using Bio;
using Bio.Algorithms.Kmer;
using RDotNet;
using RDotNet.Internals;
using HaploGrepSharp;

namespace MitoDataAssembler.Visualization
{
    public class MitochondrialAssemblyPlotMaker
    {
        //This is hard coded in R with same name
        private double outerRadius = .45;
        public int ReferenceStartCoordinate;
      
        private readonly NucMerQueryable assemblyAligner;
        private List<GenomeSlice> slices;
        private List<Arrow> arrows;
        public static Random rand = new Random();
        public readonly MitochondrialAssembly Assembly;
        /// <summary>
        /// The slices.
        /// </summary>        
		public MitochondrialAssemblyPlotMaker(MitochondrialAssembly assembly)
		{
			this.Assembly = assembly;
            slices = new List<GenomeSlice>(assembly.AllNodesInGraph.Select(x=>new GenomeSlice(x)));
            if(this.Assembly.SuccessfulAssembly)
            {
                foreach(var s in slices)
                {
                    if(Assembly.NodesInCompleteAssembly.Contains(s.Node)) {s.InMainAssembly=true;}
                }
                assemblyAligner=new NucMerQueryable(Assembly.GreedyPathAssembly.Sequence);
                assignLeftRightPositionsToSlices();
                assignHeightPositionsToSlices();
                arrows = Assembly.AllNodesInGraph.SelectMany(x => x.GetAllEdges()).Where(z => !z.IsInferiorEdge).Select(x=>new Arrow(x)).ToList();
            }
		}        
        /// <summary>
        /// Set the suggested left/right position and heights of slices based on alignment or nearby neighbors
        /// </summary>
        private void assignLeftRightPositionsToSlices()
        {
            if(slices.Count==1)
            {
                var s =slices[0];
                s.SuggestedLeftSidePosition=0;s.SuggestedRightSidePosition=Assembly.AssemblyLength;
            }
            else
            {
                //assign positions to those that can be based on nucmer, return list of remainder that can't be.
                var unAssigned=slices.Where(x=>!setPositionsBasedOnAlignment(x)).ToList();
                 //a little messy here, going to try for each unassigned slice to grab a neighbor and set values based on that
                //very inefficient at present, assuming we aren't dealing with many searches so repeated linear is fine.
                foreach(var s in unAssigned)
                {
                    var right=s.Node.GetRightEdges().Select(x=>x.ToNode.parentSlice).FirstOrDefault();
                    var left = s.Node.GetLeftEdges().Select(x => x.ToNode.parentSlice).FirstOrDefault();
                    
                    if(right!=null && right.SuggestedLeftSidePosition.HasValue) {s.SuggestedRightSidePosition=right.SuggestedLeftSidePosition.Value;}
                    if(left !=null && left.SuggestedRightSidePosition.HasValue) {s.SuggestedLeftSidePosition=left.SuggestedRightSidePosition.Value;}
                    if(!s.SuggestedLeftSidePosition.HasValue | !s.SuggestedRightSidePosition.HasValue)
                    {
                        if(s.SuggestedRightSidePosition.HasValue)
                        {
                            s.SuggestedLeftSidePosition=s.SuggestedRightSidePosition.Value-s.Node.LengthOfNode;
                        }
                        else if(s.SuggestedLeftSidePosition.HasValue) {s.SuggestedRightSidePosition=s.SuggestedLeftSidePosition+s.Node.LengthOfNode;}
                        else{s.SuggestedLeftSidePosition=GetRandomPosition();s.SuggestedRightSidePosition=s.SuggestedLeftSidePosition.Value+s.Node.LengthOfNode;}
                    }
                }               
            }
        }
        private bool setPositionsBasedOnAlignment(GenomeSlice slice)
        {
            Sequence s = new Sequence(DnaAlphabet.Instance,slice.Node.Sequence,false);
            var delts=assemblyAligner.GetDeltaAlignments(s).SelectMany(x=>x).ToList();
            if(delts.Count>0)
            {
                //check if reversed 
                if(delts[0].IsReverseQueryDirection) 
                {
                    slice.Node.ReversePath();
                    s = new Sequence(DnaAlphabet.Instance,slice.Node.Sequence,false);
                    delts=assemblyAligner.GetDeltaAlignments(s).SelectMany(x=>x).ToList();
                }
                if(delts.Count==1)
                {
                    slice.SuggestedLeftSidePosition=delts[0].FirstSequenceStart;
                    slice.SuggestedRightSidePosition=delts[0].FirstSequenceEnd;
                }
                else if(delts.Count==2) //split, happens in event of wrap around
                {
                    slice.SuggestedLeftSidePosition=delts.MaxBy(x=>x.FirstSequenceStart).FirstSequenceStart;
                    slice.SuggestedRightSidePosition=delts.MinBy(x=>x.FirstSequenceEnd).FirstSequenceEnd;
                }
                else{return false;}
                return true;
            }
            return false;
        }
        private void assignHeightPositionsToSlices()
        {
            var totalLayers=slices.Where(x=>!x.InMainAssembly).Count();
            totalLayers++;
            //now to assign heights
            double lowValue = .25;
            double highValue= outerRadius*.95;
            double delta=(highValue-lowValue)/totalLayers;
            
            double lowMain=highValue-delta;
            double height=delta*.95;
            foreach (var s in slices)
            {
                if (!s.InMainAssembly)
                {   
                    s.MagnitudeLow=lowValue;
                    s.Height=height;
                    lowValue+=delta;
                }
                else{s.MagnitudeLow=lowMain;s.Height=height;}
            }         
        }
        /// <summary>
        /// Render a PDF image to the file with the given name
        /// </summary>
        /// <param name="rInt"></param>
        /// <param name="filename"></param>
        public void Render(RInterface rInt,string filename)
        {            
//#if FALSE
            var e=rInt.CurrentEngine;
            lock(e)
            {
                if (this.slices.Count == 0 || !Assembly.SuccessfulAssembly) {
                    return;
                }
                //load the file
                string rFile = System.Reflection.Assembly.GetExecutingAssembly().Location;
				char splitChar = '\\';
				if (Bio.CrossPlatform.Environment.RunningInMono) {
					splitChar = '/';
				}
				var tempSplit=rFile.Split(splitChar);
                rFile=String.Join("/",tempSplit.Take(tempSplit.Length-1))+"/Visualization/MakeGenomeAssemblyPlot.r";
                rFile=rFile.Replace('\\','/');
                string loadFile="source('"+rFile+"')";
                e.Evaluate(loadFile);
                //TODO: THIS LOGIC REALLY NEEDS SOME CLEANING
                var curMax = Math.Max(this.slices.Max(x => x.ActualLeftGraphicalPosition.Value), this.slices.Max(x => x.ActualRightGraphicalPosition.Value));
                //e.Evaluate("assemblyLength="+Assembly.AssemblyLength.ToString());
                e.Evaluate("assemblyLength="+curMax.ToString());
                
                double outer=e.GetSymbol("outerRadius").AsNumeric().First();
                e.Evaluate("createPlot(\""+filename+"\")");
                double total = Assembly.AssemblyLength;
                //TODO: Not sure how this ever could have been possible
                if (Math.Abs(total) < double.Epsilon) {
                    return;
                }
               
                foreach (var slice in this.slices)
                {
                    Dictionary<string, object> args = new Dictionary<string, object>();
                    //drawSegment<-function(start,end,low,height)
                    args["start"] = slice.ActualLeftGraphicalPosition;
                    args["end"] = slice.ActualRightGraphicalPosition;
                    args["low"] = slice.MagnitudeLow;
                    args["height"] = slice.Height;
                    RDotNet.SymbolicExpression se;
                    
                    bool worked=e.EvaluateMethodWithNamedArguments("drawSegment", args,out se);
                    if (!worked)
                    {
                       e.Evaluate("traceback()");
                       throw new RInterfaceException("Could not draw a segment from the assembly");
                    }
                }
                //now the arrows
                //drawArrow<-function(sMag,sAng,eMag,eAng) 
                foreach (var arrow in this.arrows)
                {
                    Dictionary<string, object> args = new Dictionary<string, object>();
                    //drawSegment<-function(start,end,low,height)
                    args["sMag"] = arrow.sMag;
                    args["sAng"] = arrow.sAng;
                    args["eMag"] = arrow.eMag;
                    args["eAng"] = arrow.eAng;
                    RDotNet.SymbolicExpression se;
                    bool worked = e.EvaluateMethodWithNamedArguments("drawArrow", args, out se);
                    if (!worked)
                    {
                        throw new RInterfaceException("Could not draw an arrow from the assembly");
                    }
                }
                e.Evaluate("dev.off()");
            }
//#endif
        }



#if FALSE
               //Get a collection of all edges, and try to fill them with approximate coordinates,
                List<MetaNode.Edge> EdgesToFill = new List<MetaNode.Edge>();
                List<MetaNode.Edge> RemainingEdgeToFill = new List<MetaNode.Edge>(  
                Assembly.AllNodesInGraph.SelectMany(x=>x.GetAllEdges()));
                int PassCount= 0;
                //try to fill all positions in, after 5 passes, just assign random values
                //nodes are likely not connected to actual values.  
                while(RemainingEdgeToFill.Count>0)
                {
                    foreach(MetaNode.Edge e in RemainingEdgeToFill)
                    {
                        
                        var side1= e.GoingRight ? e.FromNode.SuggestedRightSidePosition : e.FromNode.SuggestedLeftSidePosition;
                        bool connectsOnLeft=e.GoingRight ^ e.DifferentOrientation;
                        var side2= connectsOnLeft ? e.ToNode.SuggestedLeftSidePosition : e.ToNode.SuggestedRightSidePosition;
                        //one filled, get the other too
                        if(side1.HasValue ^ side2.HasValue)
                        {
                            if(side1.HasValue)
                            {
                                if(connectsOnLeft)
                                {
                                    e.ToNode.SuggestedLeftSidePosition=side1.Value+1;
                                }
                                else{ e.ToNode.SuggestedRightSidePosition=side1.Value+1;}
                            }
                            if(side2.HasValue)
                            {
                                //for possibilities here
                                double newValue=side2.Value+1;
                                if(e.GoingRight)
                                {
                                    e.FromNode.SuggestedRightSidePosition=newValue;
                                }
                                else{ e.FromNode.SuggestedLeftSidePosition=newValue;}
                            }
                        }
                        //can we fill in from other side of node?
                        //necessary to get dead ends in, start doing on second pass 
                        else if( PassCount>1)
                        {
                            MetaNode[] nodes = new[] {e.FromNode,e.ToNode};
                            foreach(var node in nodes)
                            {
                                if(node.SuggestedRightSidePosition.HasValue && !node.SuggestedLeftSidePosition.HasValue)
                                {
                                    node.SuggestedLeftSidePosition=node.SuggestedRightSidePosition-node.LengthOfNode;
                                }
                                else if(!node.SuggestedRightSidePosition.HasValue && node.SuggestedLeftSidePosition.HasValue)
                                {
                                    node.SuggestedRightSidePosition=node.SuggestedLeftSidePosition+node.LengthOfNode;
                                }                               
                            }                           
                        }
                        else if(PassCount>5)
                        {
                            //need to randomly assign to ensure this loop exits for cliques that are disjoined from originally assigned nodes
                            if(!e.FromNode.SuggestedLeftSidePosition.HasValue)
                            {
                                e.FromNode.SuggestedLeftSidePosition = GetRandomPosition();
                            }
                        }
                        if(!side1.HasValue || !side2.HasValue)
                        {
                             EdgesToFill.Add(e);
                        }
                    }
                    PassCount++;
                    RemainingEdgeToFill=EdgesToFill;
                    EdgesToFill.Clear();
                }
                //now to fill all sides that weren't reached, e.g. those with no connections.
                //also flipped any reversed nodes
                foreach (var node in Assembly.AllNodesInGraph)
                {
                    if (!node.SuggestedLeftSidePosition.HasValue)
                    {
                        if (!node.SuggestedRightSidePosition.HasValue)
                        {
                            node.SuggestedLeftSidePosition = GetRandomPosition();
                            node.SuggestedRightSidePosition = node.SuggestedLeftSidePosition.Value + node.LengthOfNode;
                        }
                        else
                        {
                            node.SuggestedLeftSidePosition = node.SuggestedRightSidePosition.Value - node.LengthOfNode;
                        }
                    }
                    if (!node.SuggestedRightSidePosition.HasValue)
                    {
                        //should be impossible as was handled above
                        if (!node.SuggestedLeftSidePosition.HasValue)
                        {
                            node.SuggestedRightSidePosition = GetRandomPosition();
                            node.SuggestedLeftSidePosition = node.SuggestedRightSidePosition.Value - node.LengthOfNode;
                        }
                        else
                        {
                            node.SuggestedRightSidePosition = node.SuggestedLeftSidePosition.Value + node.LengthOfNode;
                        }
                    }
                    
                }
#endif

        
        private double GetRandomPosition()
        {
            return (Assembly.AssemblyLength - 2) * rand.NextDouble();
        }
       
       


       
    }
}
