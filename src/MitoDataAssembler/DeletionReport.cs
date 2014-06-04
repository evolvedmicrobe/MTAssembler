using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly;
using Bio.Algorithms.Assembly.Padena;
using Bio.Algorithms.Assembly.Graph;
using System.Diagnostics;
namespace MitoDataAssembler
{
    /// <summary>
    /// Does not allow deletions to span the control region.
    /// </summary>
    public class DeletionReport
    {
        
        private static int DeletionReportCounter=0;
        public int DeletionNumber{get;private set;}
        public bool ReferenceValuesChangeMonotonically
        {
            get;
            private set;
        }
        public List<UnAccountedForSpan> MissingSpans;
        public class UnAccountedForSpan
        {
            public static int SpanCounter=0;
            public int SubDeletionNumber{private set;get;}
            public int StartReference
            {
                get { return parent.Assembly.constituentNodes[StartNodeIndex].ReferenceGenomePosition; }
            }
            public int EndReference
            {
                get { return parent.Assembly.constituentNodes[EndNodeIndex].ReferenceGenomePosition; }
            }
            public int DeltaSequence
            {                
                get { return EndReference- StartReference- 1; }
            }
            public int UnMarkedNodesInBetween
            {  //How many nodes are missing?
                get { return EndNodeIndex - StartNodeIndex - 1; }
            }
            public int StartNodeIndex;
            public int EndNodeIndex;
            public DeletionReport parent;
            public UnAccountedForSpan(DeletionReport parent, int startNodeIndex, int endNodeIndex)
            {
                this.SubDeletionNumber=SpanCounter++;
                this.parent = parent;
                this.StartNodeIndex = startNodeIndex;
                this.EndNodeIndex = endNodeIndex;
                //now follow the path down from the start to the end see if there is one or more bifurcation.
                if (UnMarkedNodesInBetween > 0)
                {
                    int start = StartNodeIndex;
                    int next = start+1;
                    DeBruijnNode startNode = parent.Assembly.constituentNodes[start];
                    DeBruijnNode nextNode = parent.Assembly.constituentNodes[next];
                    DeBruijnNode endNode = parent.Assembly.constituentNodes[endNodeIndex];
                    bool goingLeft = startNode.GetLeftExtensionNodes().Contains(nextNode);
                    Debug.Assert(goingLeft || startNode.GetRightExtensionNodes().Contains(nextNode));
                    int foundSplitCount=0;
                    DeBruijnNode current=startNode;
                    //follow through the path finding the right guy
                    while (current!=endNode)
                    {
                        var nextNodes = goingLeft ? current.GetLeftExtensionNodesWithOrientation() : current.GetRightExtensionNodesWithOrientation();
                        if (nextNodes.Count > 1)//furcation
                        {
                            foundSplitCount++;
                            var nexti = parent.Assembly.constituentNodes[next].KmerCount;
                            decimal countNext=(decimal) parent.Assembly.constituentNodes[next].KmerCount;
                            decimal total=nextNodes.Sum(x=>x.Key.KmerCount);
                            FractionEvidence=(double) (countNext/total);
                        }
                        var nextActual = nextNodes.Where(x => x.Key == this.parent.Assembly.constituentNodes[next]).First();
                        goingLeft = !goingLeft ^ nextActual.Value;
                        next++;
                        current = nextActual.Key;
                    }
                    if (foundSplitCount == 1)
                    { SimpleBifurcation = true; }
                    else { FractionEvidence = Double.NaN; }
                }
                
            }
            public bool SimpleBifurcation = false;
            public double FractionEvidence;           
        }
        public PossibleAssembly Assembly;
        public static MitochondrialAssembly CompleteAssembly;
        
        public bool HasMutation
        {
            get {return MissingSpans.Count>1;}
        }
        public DeletionReport(PossibleAssembly assemblyToCheck )
        {
            if (CompleteAssembly == null)
            { throw new ArgumentNullException("CompleteAssembly"); }
            DeletionNumber=DeletionReportCounter++;
            MissingSpans = new List<UnAccountedForSpan>();
            Assembly = assemblyToCheck;
            LookForDeletion();
        }
        private void LookForDeletion()
        {
            bool movingUp;
            int? indexOfNodeWithLastKnownPosition=null;
            var v=Assembly.constituentNodes.Where(x=>x.IsInReference).ToList();
            var difs=Enumerable.Zip(v.Skip(1),v.Take(v.Count-1),(x,y)=> { if(x.ReferenceGenomePosition>y.ReferenceGenomePosition) return 1;
            else return 0;}).Sum();
            //only one change allowed if monotonically increasing
            if(difs>v.Count/2)
                movingUp=true;
            else
                movingUp = false;
            //now is there only one change of sign (perhaps equal to the flip?
            if(difs<2 || difs>(v.Count-1))
                ReferenceValuesChangeMonotonically=true;
            if(!movingUp)
            {
                Assembly.ReversePath();
                movingUp=true;
            }
            //Every node should be connected to a reference node, therefore all paths should be to.
            for(int i=0;i<Assembly.constituentNodes.Count;i++)
            {
                var cur=Assembly.constituentNodes[i];
                if(cur.IsInReference)
                {
                    int currentPos=cur.ReferenceGenomePosition;
                    if(indexOfNodeWithLastKnownPosition.HasValue)
                    {
                        var previous=Assembly.constituentNodes[indexOfNodeWithLastKnownPosition.Value];
                        int deltaNodes=i-indexOfNodeWithLastKnownPosition.Value;
                        int delta=cur.ReferenceGenomePosition-previous.ReferenceGenomePosition;
                        //TODO: Very ad-hoc values here
                        if(previous.ReferenceGenomePosition>15000 && currentPos<50)
                        {delta=MitoDataAssembler.Visualization.MetaNode.AssemblyLength-indexOfNodeWithLastKnownPosition.Value + currentPos;}
                        if (delta != deltaNodes)
                        {
                            UnAccountedForSpan uac = new UnAccountedForSpan(this, indexOfNodeWithLastKnownPosition.Value, i);
                            MissingSpans.Add(uac);
                        }
                    }
                    indexOfNodeWithLastKnownPosition = i;
                }
            }
        }
        public static string DeletionReportHeaderLine()
        {
            return String.Join(",", OutputColumn.OutputColumnCollection.Select(x => x.Name).ToArray());
        }
        public IEnumerable<string> DeletionReportDataLines()
        {
            //TODO Problem here
            foreach(UnAccountedForSpan span in this.MissingSpans)
                yield return String.Join(",", OutputColumn.OutputColumnCollection.Select(x => x.outFunc(span)).ToArray());
        }
    }
    public sealed class OutputColumn
    {
        public readonly string Name;
        public readonly Func<DeletionReport.UnAccountedForSpan, object> outFunc;
        public OutputColumn(string name, Func<DeletionReport.UnAccountedForSpan, object> outputFunction)
        {
            this.Name = name;
            this.outFunc = outputFunction;
        }

        public static string SafeGet(Func<DeletionReport.UnAccountedForSpan, object> func, DeletionReport.UnAccountedForSpan GC)
        {
            try
            {
                object o = func(GC).ToString();
                if (o is double)
                {
                    double n = (double)o;
                    return n.ToString("n5");
                }
                else
                {
                    return o.ToString();
                }
            }
            catch
            {
                return "NA";
            }
        }

        /// <summary>
        /// Collection of columns to output with various summary statistics.
        /// </summary>
        public static List<OutputColumn> OutputColumnCollection = new List<OutputColumn>() {
                new OutputColumn("PossibleAssemblyNumber", x => x.parent.DeletionNumber.ToString()),
                new OutputColumn("UnAccountedForSpanNumber",x=>x.SubDeletionNumber.ToString()),
                new OutputColumn("StartRef",x=>x.StartReference.ToString()),
                new OutputColumn("EndRef",x=>x.EndReference.ToString()),
                new OutputColumn("DeltaRef",x=>x.DeltaSequence.ToString()),
                new OutputColumn("StartNode",x=>x.StartNodeIndex.ToString()),
                new OutputColumn("EndNode",x=>x.EndNodeIndex.ToString()),
                new OutputColumn("SimpleBifurcation",x=>x.SimpleBifurcation.ToString()),
                new OutputColumn("FractionEvidence",x=>x.FractionEvidence.ToString()),
                new OutputColumn("RefChangesMonotonically",x=>x.parent.ReferenceValuesChangeMonotonically.ToString()),
                new OutputColumn("AssemblyLength",x=>x.parent.Assembly.DirtySequence.Count.ToString()),
                new OutputColumn("SequenceOfAssembly",x=>x.parent.Assembly.DirtySequence.ConvertToString())

        };
    }
}
