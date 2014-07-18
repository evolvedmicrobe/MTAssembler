using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly;
using Bio.Algorithms.Assembly.Graph;
using Bio;
using Bio.Algorithms.Kmer;
using System.Diagnostics;
using MitoDataAssembler.Extensions;
using MitoDataAssembler.Utils;

namespace MitoDataAssembler.Visualization
{      
    /// <summary>
    /// A Node that represents a contiguous path in a graph
    /// (that is several subnodes) that occur with no divergence. Used for building visualizations and assemblies of the graph.
    /// </summary>
    [DebuggerDisplay("Value = {NodeNumber}")]
    public class MetaNode
    {
         /// <summary>
        /// The K-mer length used for the data
        /// </summary>
        public static int KmerLength;
        //used to get information about the parent, somewhat 
        //TODO: Fix this somewhat strange hack.
        public GenomeSlice parentSlice;
        public class Edge
        {
            public MetaNode ToNode; public MetaNode FromNode;
            public uint Weight; public bool DifferentOrientation;
            public bool GoingRight;
            /// <summary>
            /// Every path from A->B can also be expressed as B->A,
            /// this checks the two debruijin nodes that form the connection, and determines if the 
            /// FromNode has a greater value than the ToNode, returning true or false.
            /// </summary> 
            public bool IsInferiorEdge
            {
                get {
                    return FromInnerNode.NodeValue.KmerData > ToInnerNode.NodeValue.KmerData;
                }
            }
            public DeBruijnNode FromInnerNode
            {
                get
                {
                    if (GoingRight)
                        return FromNode.ConstituentNodes.Last();
                    else
                        return FromNode.ConstituentNodes[0];
                }
            }
            public DeBruijnNode ToInnerNode
            {
                get
                {
                    bool grabLeft = GoingRight ^ DifferentOrientation;
                    if (grabLeft)
                        return ToNode.ConstituentNodes[0];
                    else
                        return ToNode.ConstituentNodes.Last();                
                }
            }
        }
        public class PossiblePath
        {
            public MetaNode NeighborNode;
            public uint Weight;
            public bool NeedsReversing;
        }
        /// <summary>
        /// Gets or sets the default fill color.
        /// </summary>
        /// <value>The default fill color.</value>
        

        public bool CircularLoop = false;
       
        /// <summary>
        /// Classifies the node for purposes of creating a graph where links of singly joined nodes are combined
        /// </summary>
        /// <returns></returns>
        public static NODE_TYPE ClassifyNode(DeBruijnNode startNode)
        {
            var lefts = startNode.GetLeftExtensionNodes().ToArray();
            var rights = startNode.GetRightExtensionNodes().ToArray();
            if (lefts.Any(x => rights.Contains(x)))
            {
                return NODE_TYPE.END_LOOPS_ON_ITSELF;
            }
            //First to check if this guy can form an infinite circle with itself
            int validLeftExtensionsCount = lefts.Length; 
            int validRightExtensionsCount = rights.Length;
            if (validLeftExtensionsCount != 1 && validRightExtensionsCount == 1)
            {
                return NODE_TYPE.GO_RIGHT;
            }
            else if (validLeftExtensionsCount == 1 && validRightExtensionsCount != 1)
            {
                return NODE_TYPE.GO_LEFT;
            }
            else if (validRightExtensionsCount == 1 && validLeftExtensionsCount == 1)
            {   
                return NODE_TYPE.LINK_IN_CHAIN;
            }
            else if (validRightExtensionsCount > 1 && validLeftExtensionsCount > 1)
            {
                return NODE_TYPE.NEXUS;
            }
            else if (validLeftExtensionsCount != 1 && validRightExtensionsCount != 1)
            {
                return NODE_TYPE.ISLAND;
            }
            throw new Exception("Apparently you did not handle all cases...");
        }
        
        public enum NODE_TYPE
        {
            END_LOOPS_ON_ITSELF,//A node can be traced to refer to itself infitinely on one side A->A->A->A
            GO_LEFT, // Right has no or >1 edge
            GO_RIGHT, // Left has no or >1 edge
            ISLAND, // No edges on either side
            LINK_IN_CHAIN,//single nodes on either side
            NEXUS,//silly name for node with >1 on each side
        };
        //Both sides have 1
        public List<DeBruijnNode> ConstituentNodes=new List<DeBruijnNode>();
        //public List<Edge> RightNodes = new List<Edge>();
        //public List<Edge> LeftNodes = new List<Edge>();
        List<byte> contigSequence;
        public double AvgKmerCoverage;       
       
        public int LengthOfNode
        {
            get { return this.Sequence.Length;}
        }

        public string ReverseComplementedSequence
        {
            get
            {
                var tmpSequence = new Sequence(DnaAlphabet.Instance, Sequence);
                tmpSequence = new Sequence(tmpSequence.GetReverseComplementedSequence());
                return tmpSequence.ConvertToString(0, tmpSequence.Count);
            }
        }

       
        public void ReversePath()
        {
            ConstituentNodes.Reverse();
            var tmpSequence = new Sequence(DnaAlphabet.Instance, contigSequence.ToArray());
            tmpSequence = new Sequence(tmpSequence.GetReverseComplementedSequence());
            contigSequence = tmpSequence.ToList();
            Sequence = tmpSequence.ConvertToString(0, tmpSequence.Count);
        }
        /// <summary>
        /// Takes a node that connects to this meganode at either the top or bottom, and gives the nodes on the other side that lead away
        /// </summary>
        /// <param name="incomingNode"></param>
        /// <returns></returns>
        public IEnumerable<PossiblePath> GetOutgoingNodes()
        {
            foreach (Edge e in this.GetRightEdges())
            {
                yield return new PossiblePath()
                {
                    NeighborNode = e.ToNode as MetaNode,
                    NeedsReversing = e.DifferentOrientation,
                    Weight=e.Weight
                };
            }            
        }

        public string TrailingKmer
        {
            get
            {
                return this.Sequence.Substring(Sequence.Length-KmerLength+1,KmerLength-1);
            }
        }
        public string LeadingKmer
        {
            get
            {
                return this.Sequence.Substring(0,KmerLength-1);
            }
        }
      
        public static string ReverseComplimentSequence(string Sequence)
        {
            return (new Sequence(DnaAlphabet.Instance, Sequence).GetReverseComplementedSequence() as Sequence).ConvertToString();
        }
        public IEnumerable<Edge> GetRightEdges()
        {
            var SearchSequence = this.TrailingKmer;
            var last = ConstituentNodes.Last();
            foreach (var neighbor in GetNodesLeavingBottom())// nextNodes)
            {
                var neighborMetaNode = neighbor.ParentMetaNode as MetaNode;
                //TODO: The sequence can be palindromic, allowing the same k-mer to appear in the front and back.
                bool startsFront=neighborMetaNode.Sequence.StartsWith(SearchSequence);
                bool startsBack = neighborMetaNode.ReverseComplementedSequence.StartsWith(SearchSequence);
                    if (startsFront && !startsBack)
                    {
                        yield return new Edge() { FromNode=this, ToNode = neighborMetaNode, Weight = CalculateConnectionWeight(last, neighborMetaNode.ConstituentNodes[0]), DifferentOrientation = false,GoingRight=true };
                    }
                    else if(startsBack && !startsFront)
                    {
                        yield return new Edge() { FromNode=this, ToNode = neighborMetaNode, Weight = CalculateConnectionWeight(last, neighborMetaNode.ConstituentNodes.Last()), DifferentOrientation = true, GoingRight=true };
                    }
                    else
                    {
                              yield return new Edge() { FromNode=this, ToNode = neighborMetaNode, Weight = CalculateConnectionWeight(last, neighborMetaNode.ConstituentNodes[0]), DifferentOrientation = false,GoingRight=true };
                              yield return new Edge() { FromNode=this, ToNode = neighborMetaNode, Weight = CalculateConnectionWeight(last, neighborMetaNode.ConstituentNodes.Last()), DifferentOrientation = true, GoingRight=true };
#if !DEBUG
    throw new Exception("GRR!!! a meganode connects to both the front and back of another meganode, palindrome problems.");
#endif
                    }
             }
        }
        public IEnumerable<Edge> GetLeftEdges()
        {
            var first = this.ConstituentNodes[0];
            var SearchSequence = this.LeadingKmer;
            foreach (var neighbor in GetNodesLeavingTop())
            {
                var neighborMegaNode = neighbor.ParentMetaNode as MetaNode;
                if (neighborMegaNode.Sequence.EndsWith(SearchSequence))
                {
                    yield return new Edge()
                    {
                        GoingRight = false,
                        FromNode = this,
                        ToNode = neighborMegaNode,
                        Weight = CalculateConnectionWeight(first, neighborMegaNode.ConstituentNodes.Last()),
                        DifferentOrientation = false
                    };
                }
                else
                {
                    Debug.Assert(neighborMegaNode.ReverseComplementedSequence.EndsWith(SearchSequence));
                    yield return new Edge()
                    {
                        GoingRight = false,
                        FromNode = this,
                        ToNode = neighborMegaNode,
                        Weight = CalculateConnectionWeight(first, neighborMegaNode.ConstituentNodes.First()),
                        DifferentOrientation = true
                    };
                }
            }
        }
        public IEnumerable<DeBruijnNode> GetNodesLeavingBottom()
        {
                //if chain is longer than one, use previous node to get latest
                if (this.ConstituentNodes.Count > 1)
                {
                    var bottomNode = ConstituentNodes.Last();
                    DeBruijnNode penUltimate = ConstituentNodes[ConstituentNodes.Count - 2];
                    bool goingRight = penUltimate.GetRightExtensionNodes().Contains(bottomNode);
                    var next = goingRight ? penUltimate.GetRightExtensionNodesWithOrientation().Where(x => x.Key == bottomNode).First() :
                    penUltimate.GetLeftExtensionNodesWithOrientation().Where(x => x.Key == bottomNode).First();
                    var nextSet = goingRight ^ next.Value ? next.Key.GetLeftExtensionNodes() :
                    next.Key.GetRightExtensionNodes();
                    foreach (var k in nextSet)
                    {
                        yield return k;
                    }
                }
                else
                {
                    var baseNode = this.ConstituentNodes[0];
                    Debug.Assert(KmerLength == Sequence.Length);
                    var ns = new Sequence(DnaAlphabet.Instance, baseNode.GetOriginalSymbols(MetaNode.KmerLength));
                    bool orientationRight;// = baseNode.GetOriginalSymbols(KmerLength).SequenceEqual(new DnaAlphabet(DnaAlphabet.Instance, Sequence));
                    if (ns.ConvertToString().Equals(Sequence))
                        orientationRight = true;
                    else if ((new Sequence(ns.GetReverseComplementedSequence()).ConvertToString().Equals(Sequence)))
                        orientationRight = false;
                    else
                        throw new Exception("AAA");
                    var nextNodes = orientationRight ? baseNode.GetRightExtensionNodes() : baseNode.GetLeftExtensionNodes();
                    foreach (var v in nextNodes)
                    {
                        yield return v;
                    }
                }
            
        }
        public IEnumerable<DeBruijnNode> GetNodesLeavingTop()
        {
                //if chain is longer than one, use previous node to get latest
                if (this.ConstituentNodes.Count > 1)
                {
                    var topNode = ConstituentNodes[0];
                    DeBruijnNode penUltimate = ConstituentNodes[1];
                    bool goingLeft = penUltimate.GetLeftExtensionNodes().Contains(topNode);
                    var next = goingLeft ? penUltimate.GetLeftExtensionNodesWithOrientation().Where(x => x.Key == topNode).First() :
                    penUltimate.GetRightExtensionNodesWithOrientation().Where(x => x.Key == topNode).First();
                    var nextSet = goingLeft ^ next.Value ? next.Key.GetRightExtensionNodes() :
                    next.Key.GetLeftExtensionNodes();
                    foreach (var k in nextSet)
                    {
                        yield return k;
                    }
                }
                else
                {
                    var baseNode = this.ConstituentNodes[0];
                    Debug.Assert(KmerLength == Sequence.Length);
                    var ns = new Sequence(DnaAlphabet.Instance, baseNode.GetOriginalSymbols(MetaNode.KmerLength));
                    bool orientationRight;// = baseNode.GetOriginalSymbols(KmerLength).SequenceEqual(new DnaAlphabet(DnaAlphabet.Instance, Sequence));
                    if (ns.ConvertToString().Equals(Sequence))
                        orientationRight = true;
                    else if ((new Sequence(ns.GetReverseComplementedSequence()).ConvertToString().Equals(Sequence)))
                        orientationRight = false;
                    else
                        throw new Exception("AAA");
                    var nextNodes = orientationRight ? baseNode.GetLeftExtensionNodes() : baseNode.GetRightExtensionNodes();
                    foreach (var v in nextNodes)
                    {
                        yield return v;
                    }
                }            
        }
        public IEnumerable<Edge> GetAllEdges()
        {
            foreach (Edge e in GetLeftEdges())
            {
                yield return e;
            }
            foreach (Edge e in GetRightEdges())
            {
                yield return e;
            }
        }
        
        public string Sequence;
        public int NodeNumber;
        
       
        public static uint CalculateConnectionWeight(DeBruijnNode FirstNode, DeBruijnNode SecondNode)
        {
            //First verify that they share
            if (!FirstNode.GetExtensionNodes().Contains(SecondNode))
            {
                throw new Exception("Can't calculate non-overlapping extensions");
            }
            return SecondNode.KmerCount;
        }
        public MetaNode() { }
        public MetaNode(DeBruijnNode startNode, DeBruijnGraph graph)
        {
            this.NodeNumber = GraphGenerator.NodeCount++;
            KmerLength=graph.KmerLength;            
            if (startNode.IsVisited)
            {throw new Exception("If a node has been visited it should not form a metanode, suggests an infinite recursion problem");}
            NODE_TYPE type = ClassifyNode(startNode);
            startNode.IsVisited=true;
            //Either of these become their own thing
            if (type == NODE_TYPE.NEXUS || type == NODE_TYPE.ISLAND || type==NODE_TYPE.END_LOOPS_ON_ITSELF)
            {
                ConstituentNodes.Add(startNode);
                contigSequence = new List<byte>(graph.GetNodeSequence(startNode));
                Sequence = (new Sequence((IAlphabet)NoGapDnaAlphabet.Instance, contigSequence.ToArray())).ConvertToString(0, contigSequence.Count);
            }
            else if (type == NODE_TYPE.LINK_IN_CHAIN)
            {
                contigSequence = new List<byte>(graph.GetNodeSequence(startNode));
                if (!VerifyNotCircular(startNode))
                {
                    MakeCircle(startNode, graph);
                    //throw new Exception("Non circular visualizations not currently supported");
                }
                else
                {
                    //go right first
                    contigSequence = new List<byte>(graph.GetNodeSequence(startNode));
                    //var nextNodes = ExtendChain(startNode, true, graph);
                    ExtendChain(startNode, true, graph);
                    //copy the right information and clear it out
                    var tmpRightSeq = contigSequence.ToArray();
                    //skip the first node
                    var tmpRightNodes = ConstituentNodes.Skip(1).ToArray();
                    ConstituentNodes.Clear();
                    contigSequence.Clear();
                    //now go left
                    ExtendChain(startNode, false, graph);
                    //now lets combine
                    ConstituentNodes.Reverse();
                    ConstituentNodes.AddRange(tmpRightNodes);
                    var tmpSequence = new Sequence(DnaAlphabet.Instance, contigSequence.ToArray());
                    tmpSequence = new Sequence(tmpSequence.GetReverseComplementedSequence());
                    string LeftSequence = "";
                    if (tmpSequence.Count > 0)
                    {
                        LeftSequence = tmpSequence.ConvertToString(0, tmpSequence.Count);
                    }
                    tmpSequence = new Sequence(DnaAlphabet.Instance, tmpRightSeq);
                    Sequence = LeftSequence + tmpSequence.ConvertToString(0, (tmpSequence.Count));
                    contigSequence = new Sequence(DnaAlphabet.Instance, Sequence).ToList();
                }
            }
            else if (type == NODE_TYPE.GO_LEFT)
            {
                contigSequence = new List<byte>(graph.GetNodeSequence(startNode).GetReverseComplementedSequence());
                //var nextNodes = ExtendChain(startNode, false, graph);
                ExtendChain(startNode, false, graph);
                var tmpSequence = new Sequence(DnaAlphabet.Instance, contigSequence.ToArray());
                //somewhat confusing - originally built the RC of sequence, so RCing again to get correct orientation for 
                //neighbors
                
                tmpSequence = new Sequence(tmpSequence.GetReverseComplementedSequence());
                contigSequence = tmpSequence.ToList();
                Sequence = tmpSequence.ConvertToString(0, tmpSequence.Count);
                //flip it so nodes and sequence are in order
                ConstituentNodes.Reverse();
                
            }
            else if (type == NODE_TYPE.GO_RIGHT)
            {
                contigSequence = new List<byte>(graph.GetNodeSequence(startNode));
                //var nextNodes = ExtendChain(startNode, true, graph);
                ExtendChain(startNode, true, graph);
                var tmpSequence = new Sequence(DnaAlphabet.Instance, contigSequence.ToArray());
                Sequence = tmpSequence.ConvertToString(0, tmpSequence.Count);
            }
            
            Cement();

        }

        /// <summary>
        /// This gets the next symbol from a node while forming chains.  This can be made a lot more efficient if it turns in to a bottleneck.
        /// all chains are extended from either the first or last base present in the node, and this base is either forward
        /// or reverse complimented, this method reflects this.
        /// </summary>
        /// <param name="node">Next node</param>
        /// <param name="graph">Graph to get symbol from</param>
        /// <param name="GetFirstNotLast">First or last base?</param>
        /// <param name="ReverseComplimentBase">Should the compliment of the base be returned</param>
        /// <returns></returns>
        public static byte GetNextSymbol(DeBruijnNode node,int kmerLength, bool GetRCofFirstBaseInsteadOfLastBase)
        {
            if (node == null)
            {throw new ArgumentNullException("node");}
            byte[] symbols=node.GetOriginalSymbols(kmerLength);
            byte value = GetRCofFirstBaseInsteadOfLastBase ? symbols.First() : symbols.Last();
            if (GetRCofFirstBaseInsteadOfLastBase)
            {
                byte value2;
                bool rced=DnaAlphabet.Instance.TryGetComplementSymbol(value,out value2);
                //Should never happend
                if(!rced)
                {
                    throw new Exception("Could not revcomp base during graph construction");
                }
                value=value2;
            }            
            return value;
           
        }
       

        /// <summary>
        /// Follow a chain along a path link a bifurcation or no additional nodes appear.   
        /// </summary>
        /// <param name="currentNode"></param>
        /// <param name="goRight"></param>
        /// <param name="graph"></param>
        /// <returns></returns>
        private Dictionary<DeBruijnNode,bool> ExtendChain(DeBruijnNode currentNode, bool goRight,DeBruijnGraph graph)
        {
            ConstituentNodes.Add(currentNode);
            currentNode.IsVisited=true;
            Dictionary<DeBruijnNode, bool> nextNodes;
            if (goRight)
            {nextNodes = currentNode.GetRightExtensionNodesWithOrientation();}
            else
            {nextNodes = currentNode.GetLeftExtensionNodesWithOrientation();}
            DeBruijnNode next;
            DeBruijnNode last=currentNode;
            while (nextNodes.Count == 1)
            {
                var nextSet = nextNodes.First();
                next = nextSet.Key;
                bool sameOrientation = nextSet.Value;
                goRight = (!goRight) ^ sameOrientation;
                int oppositeDirectionExtensions = goRight ? next.LeftExtensionNodesCount : next.RightExtensionNodesCount;
                int sameDirectionExtensions = goRight ? next.RightExtensionNodesCount : next.LeftExtensionNodesCount;
                Debug.Assert(oppositeDirectionExtensions != 0);//should always be >1 given the node we came from.
                if (oppositeDirectionExtensions > 1)
                {
                    break;//nexus, or need to start a new node, no visit count
                }
                else
                {
                    //we have to check if the right path loops back on itself, for example TTTTTTTTTTTT could keep adding T's to infinity, always going back to the same node.
                    //However, it is also possible that the node can refer to itself, but not in a loop, e.g. by turning around, like
                    //TTTTTTTCAATTGAAAAAA which matches the reverse compliment of itself, so leaves the other side (not this might be incorrect as this is guaranteed).
                    //unfortunately, impossible to tell without looking two steps into the future, and because we are doing this one at a time,
                    //have to unwind the last addition.
                    if (next.IsVisited)
                    {
                        //note that this is a bit of an unusual step, as most of the time the other direction extensions will be >1.  This can only 
                        //happen if the only incoming node to this k-mer-1 palindrome does not have any other links, which will be rare.
                        if (next == last)
                        {
                            //if going to refer to itself again, it's a loop, need to end it and make a new self referencing mega node.
                            var temp = goRight ? next.GetRightExtensionNodesWithOrientation() : next.GetLeftExtensionNodesWithOrientation();
                            if (temp.Count == 1 && temp.First().Key == last)//three times in a row, need to remove this node from the list as we are not leaving in a different direction, //and need to unvisit the node
                            {
                                //unwind the last addition, this node needs to be a self-referencing mega node
                                next.IsVisited = false;
                                Debug.Assert(ConstituentNodes.Last() == next);
                                ConstituentNodes.RemoveAt(ConstituentNodes.Count - 1);
                                contigSequence.RemoveAt(ConstituentNodes.Count - 1);
                                Debug.Assert(ConstituentNodes.Last() != next);

                                //exit, we are as low as we can go.
                                break;
                            }
                            //criteria is that the sequence can't be there more than once

                        }
                        //At most a kmer can be used to represent the forward and reverse sequence that it has.
                        Debug.Assert(this.ConstituentNodes.Count(x => x == next) < 3);
                    }
                    byte nextSymbol = GetNextSymbol(next, graph.KmerLength, !goRight);
                    contigSequence.Add(nextSymbol);
                    //byte[] original=next.NodeValue.GetOriginalSymbols(MegaNode.KmerLength);
                    //var s=new Sequence(DnaAlphabet.Instance,original);
                    //Console.WriteLine(s.ConvertToString());

                    next.IsVisited = true;
                    ConstituentNodes.Add(next);
                    nextNodes = goRight ? next.GetRightExtensionNodesWithOrientation() : next.GetLeftExtensionNodesWithOrientation();
                    last = next;
                }
            }                
            
            return nextNodes;
        }      
        /// <summary>
        /// Follow a node with one neighbor on either side and make sure it never reaches itself, which is problematic for making these things.   
        /// Note that nodes can go to A->A->C if they refer to themselves but match the reverse compliment of themselves
        /// </summary>
        /// <param name="currentNode"></param>
        /// <param name="goRight"></param>
        /// <param name="graph"></param>
        /// <returns></returns>
        private bool VerifyNotCircular(DeBruijnNode currentNode)
        {
            List<DeBruijnNode> visitedNodes=new List<DeBruijnNode>();

            if(ClassifyNode(currentNode)!=NODE_TYPE.LINK_IN_CHAIN)
            {
                throw new Exception("Node type doesn't match well!");
            }
            else
            {
                //go right, if we wind up where we started, circle.  
                var nextNode=currentNode.GetRightExtensionNodesWithOrientation().First();
                bool goingRight=true;
                //we now either have the second or third node in path as next
                while(ClassifyNode(nextNode.Key)==NODE_TYPE.LINK_IN_CHAIN)
                {
                    visitedNodes.Add(nextNode.Key);
                    //determine if this is a kink or not, which will trigger issue at only first node.
                    if (nextNode.Key == currentNode)
                    {
                        //only one way to get back to the start, either we are in a circle, or the first node loops in to its reverse compliment and exits
                        //the other way, a "kink" so to speak, we know we have visited the right node since we started there, if we visited the left, problems
                        bool leftVisited = visitedNodes.Contains(currentNode.GetLeftExtensionNodes().First());
                        if (leftVisited)
                        {
                            return false;
                        }
                        Debug.Assert(visitedNodes.Contains(currentNode.GetRightExtensionNodes().First()));
                    }
                        
                    goingRight = !(goingRight ^ nextNode.Value);
                    var nextSet = goingRight ? nextNode.Key.GetRightExtensionNodesWithOrientation() : nextNode.Key.GetLeftExtensionNodesWithOrientation();
                    if (nextSet.Count != 1)
                    {
                        return true;
                    }
                    nextNode = nextSet.First();
                    
                }
                return true;
            }
        }
        private void MakeCircle(DeBruijnNode startNode,DeBruijnGraph graph)
        {
            CircularLoop = true;
            byte[] v= startNode.GetOriginalSymbols(graph.KmerLength);
            Console.WriteLine((new Sequence(DnaAlphabet.Instance, v)).ToString());
            ConstituentNodes.Add(startNode);
            startNode.IsVisited = true;
            Dictionary<DeBruijnNode, bool> nextNodes;
            bool goRight = true;
            nextNodes = startNode.GetRightExtensionNodesWithOrientation();
            var nextSet = nextNodes.First();
            DeBruijnNode next = nextSet.Key;
            while (next!=startNode)
            {
                next.IsVisited = true;
                ConstituentNodes.Add(next);
                bool sameOrientation = nextSet.Value;
                NODE_TYPE nextType = ClassifyNode(next);
                //what direction do we get the node following the next one from? (Note path out determined by path in, so don't need to look at next node to get side of the one after).
                goRight = (!goRight) ^ sameOrientation;
                if (nextType == NODE_TYPE.LINK_IN_CHAIN)
                {
                    //NOTE: four possibilities condense in to 2 possible sides so written with ^ operator
                    nextNodes = goRight ? next.GetRightExtensionNodesWithOrientation() : next.GetLeftExtensionNodesWithOrientation();
                    //now how to determine what base to get? This only depends on relationship of current node to next node
                    //in all cases we either grab the RC of the first base or the last base, and which to grab is determined by incoming node
                    byte nextSymbol = GetNextSymbol(next, graph.KmerLength, !goRight);
                    contigSequence.Add(nextSymbol);
                }
                else
                { throw new Exception("Non circular path being treated like one"); }
                nextSet = nextNodes.First();
                next = nextSet.Key;               
            }
            Sequence = (new Sequence((IAlphabet)NoGapDnaAlphabet.Instance, contigSequence.ToArray())).ConvertToString(0, contigSequence.Count);
            
        }
        public int Lowest_Reference_Position;
        public int Highest_Reference_Position;
        
        /// <summary>
        /// Rotate the node so that its top values are lower on the reference genome than values further down.  
        /// Note: Goofy name, meant to indicate we are done making the class, probably should be a factory/constructor class.
        /// </summary>
        private void Cement()
        {
            Debug.Assert(this.ConstituentNodes.Count != 0);
            Debug.Assert(this.contigSequence.Count != 0);
            this.ConstituentNodes[0].ParentMetaNode = this;
            this.ConstituentNodes.Last().ParentMetaNode = this;
            if (ConstituentNodes.Count > 0 && ConstituentNodes.Any(x => x.IsInReference))
            {
                //first to flip values if need be.
                short low = (short)ConstituentNodes.Where(x => x.IsInReference).First().ReferenceGenomePosition;
                short high = (short)ConstituentNodes.Where(x => x.IsInReference).Last().ReferenceGenomePosition;
                short ahigh = Math.Max(high, low);
                short alow = Math.Min(high, low);
                bool wrapsAround = Utils.CircularGenomeCaseHandlers.SegmentLoopsAround(alow, ahigh, LengthOfNode);
                bool needFlip = (high < low && !wrapsAround) || wrapsAround && low < high;
                if (needFlip)
                {
                    short tmp = low;
                    low = high;
                    high = tmp;
                    this.ReversePath();
                }
                Lowest_Reference_Position = low;
                Highest_Reference_Position = high;
                int lastVal = low;
                //verify a monotonic increase in sequences
                if (!this.CircularLoop)
                {
                    for (int i = 0; i < (ConstituentNodes.Count - 1); i++)
                    {
                        if (ConstituentNodes[i].IsInReference)
                        {
                            if (ConstituentNodes[i].ReferenceGenomePosition < lastVal)
                            {
                                if (wrapsAround && lastVal > 15000 && ConstituentNodes[i].ReferenceGenomePosition < 100)
                                {
                                    lastVal = ConstituentNodes[i].ReferenceGenomePosition;
                                    continue;
                                }
                                else
                                {
                                    var vals = ConstituentNodes.Select(x => x.ReferenceGenomePosition);
                                    throw new Exception("The constituent nodes do not monotonically increase along the genome, suggests a problem" + vals.Count().ToString());
                                }
                            }
                            lastVal = ConstituentNodes[i].ReferenceGenomePosition;
                        }
                    }
                }
            }
            this.AvgKmerCoverage = ConstituentNodes.Average(x => (double)x.KmerCount);
        }
    }
       
}