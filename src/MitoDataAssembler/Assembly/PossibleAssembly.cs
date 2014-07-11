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
using HaploGrepSharp;
using MitoDataAssembler.Extensions;
using MitoDataAssembler.Visualization;

namespace MitoDataAssembler
{
    /// <summary>
    /// A Node that represents a contiguous path in a graph
    /// (that is several subnodes) that occur with no divergence. Used for building visualizations and assemblies of the graph.
    /// </summary>
    public class PossibleAssembly : IList<DeBruijnNode>
    {
        private bool finalized=false;
        public bool CircularLoop = false;
        //TODO: Remove these
        public int? FirstReferencePosition
        {
            get
            {
                var v= constituentNodes.Where(x=>x.IsInReference).Select(x=>x.ReferenceGenomePosition).ToList();
                if(v.Count>0)
                    return v.First();
                else
                    return null;
            }
        }
        public int? LastReferencePosition
        {
            get
            {
                var v= constituentNodes.Where(x=>x.IsInReference).Select(x=>x.ReferenceGenomePosition).ToList();
                if(v.Count>0)
                    return v.Last();
                else
                    return null;
            }
        }
        
        private List<DeBruijnNode> constituentNodes;
        public IReadOnlyList<DeBruijnNode> ConstitutentNodes
        {
            get { return constituentNodes.AsReadOnly(); }
        }

        
        private HashSet<DeBruijnNode> nodeHash;
        protected List<byte> contigSequence;

        const int NOTSETFLAG = -999;
        private int startGenomePosition = NOTSETFLAG;

        public PossibleAssembly(DeBruijnNode start, bool startGoingRight)
        {
            this.contigSequence = new List<byte>(); 
            this.constituentNodes = new List<DeBruijnNode>();
            this.constituentNodes.Add(start);
            this.nodeHash = new HashSet<DeBruijnNode>();
            if(startGoingRight)
            {
                contigSequence = new List<byte>(LargeDeletionFinder.graph.GetNodeSequence(start));
            }
            else
            {
                contigSequence = new List<byte>(LargeDeletionFinder.graph.GetNodeSequence(start).GetReverseComplementedSequence());
            }            
        }

        /// <summary>
        /// Combines two assemblies derived from going different directions from the same start node.
        /// </summary>
        /// <param name="leftSide"></param>
        /// <param name="rightSide"></param>
        public PossibleAssembly(PossibleAssembly leftSide, PossibleAssembly rightSide)
        {
             //First to verify that both start in same location.
             if(leftSide.constituentNodes[0]!=rightSide.constituentNodes[0])
             {
                 throw new Exception("Cannot combine assemblies that start in different locations!");
             }
             var tmpRightSeq = rightSide.contigSequence.ToArray();
             //skip the first node
             var tmpRightNodes = rightSide.constituentNodes.Skip(1).ToArray();
             constituentNodes=leftSide.constituentNodes.ToList();
             constituentNodes.Reverse();
             //now lets combine
             constituentNodes.AddRange(tmpRightNodes);
             var tmpSequence = new Sequence(DnaAlphabet.Instance, leftSide.contigSequence.ToArray());
             tmpSequence = new Sequence(tmpSequence.GetReverseComplementedSequence());
             string LeftSequence = "";
             LeftSequence = tmpSequence.ConvertToString(0, tmpSequence.Count);
             tmpSequence = new Sequence(DnaAlphabet.Instance, tmpRightSeq);
             string tmpSequence2 = LeftSequence + tmpSequence.ConvertToString(LargeDeletionFinder.graph.KmerLength, (tmpSequence.Count-LargeDeletionFinder.graph.KmerLength));
             contigSequence = new Sequence(DnaAlphabet.Instance, tmpSequence2).ToList();
         }

        public double AvgKmerCoverage
         {
             get {
                 return constituentNodes.Average(x=>(double)x.KmerCount);
             }
         }
        public void ReversePath()
        {
            constituentNodes.Reverse();
            var tmpSequence = new Sequence(DnaAlphabet.Instance, contigSequence.ToArray());
            tmpSequence = new Sequence(tmpSequence.GetReverseComplementedSequence());
            contigSequence = tmpSequence.ToList();
            
        }
        public string SequenceAsString
        {
            get
            {
                return this.Sequence.ConvertToString();
            }
        }
        /// <summary>
        /// Gets the current sequence, not a finalized version that is correctly oriented
        /// </summary>
        public Sequence DirtySequence
        {
            get {
                if (finalized) return _sequence;
                else
                    return new Sequence(DnaAlphabet.Instance, this.contigSequence.ToArray(), false); }
        }
        private Sequence _sequence;
        public Sequence Sequence
        {
            get
            {

                if (!finalized) { throw new Exception("Can't get sequence before one is finished making the possible assembly"); }
                else
                {
                    return _sequence;
                }
            }
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
        private byte GetNextSymbol(DeBruijnNode node, DeBruijnGraph graph, bool GetRCofFirstBaseInsteadOfLastBase)
        {
            if (node == null)
            { throw new ArgumentNullException("node"); }
            int kmerLength = graph.KmerLength;
            byte[] symbols = node.GetOriginalSymbols(graph.KmerLength);
            byte value = GetRCofFirstBaseInsteadOfLastBase ? symbols.First() : symbols.Last();
            if (GetRCofFirstBaseInsteadOfLastBase)
            {
                byte value2;
                bool rced = DnaAlphabet.Instance.TryGetComplementSymbol(value, out value2);
                //Should never happend
                if (!rced)
                {
                    throw new Exception("Could not revcomp base during graph construction");
                }
                value = value2;
            }
            return value;
        }

          
        /// <summary>
        /// Adds the node to the path, if the node is a painted node, returns a value 
        /// indicating the distance between the two
        /// </summary>
        /// <param name="newNode"></param>
        /// <returns></returns>
        public void Add(DeBruijnNode newNode, byte symbolFromNode)
        {
            constituentNodes.Add(newNode);
            nodeHash.Add(newNode);
            contigSequence.Add(symbolFromNode);
        }
      
        public PossibleAssembly()
        {
            constituentNodes=new List<DeBruijnNode>();
            contigSequence = new List<byte>();
            nodeHash = new HashSet<DeBruijnNode>();
        }
        public PossibleAssembly(List<DeBruijnNode> nodes, List<byte> contigSequence)
        {

            this.constituentNodes = nodes.ToList();
            this.contigSequence = contigSequence.ToList();
            nodeHash = new HashSet<DeBruijnNode>(nodes);
        }
        /// <summary>
        /// Create a deep copy of the list and sorts it so that it is easy to identify redundant paths node comes first
        /// </summary>
        /// <returns></returns>
        public PossibleAssembly Clone()
        {
            PossibleAssembly pdb = new PossibleAssembly(this.constituentNodes,this.contigSequence);
            return pdb;
        }

        /// <summary>
        /// Adds a meganode to the current assembly by adding it's constituent nodes and
        /// discarding any information form the meganode itself.
        /// </summary>
        /// <param name="nodeToAdd"></param>
        public void AddMetaNode(MetaNode nodeToAdd)
        {
            foreach (var node in nodeToAdd.ConstituentNodes)
            {
                this.constituentNodes.Add(node);
            }
            int kmerLength=nodeToAdd.LeadingKmer.Length;
            this.contigSequence.AddRange(new Sequence(DnaAlphabet.Instance,nodeToAdd.Sequence.Substring(kmerLength, nodeToAdd.Sequence.Length - kmerLength))); 
        }
        /// <summary>
        /// Orient this sequence to the reference, give its MSA with it. 
        /// </summary>
        public void FinalizeAndOrientToReference()
        {
            if (!finalized)
            {
                _sequence = new Bio.Sequence(DnaAlphabet.Instance, contigSequence.ToArray());
                var delts = ReferenceGenome.GetDeltaAlignments(_sequence).SelectMany(x => x).ToList();
                //now to change orientation so that it matches well with the reference
                //will do a simple single flip
                if (delts.Count > 1)
                {
                    // Find the lowest value on the reference, go from there
                    var deltToUse = delts.MinBy(x => x.FirstSequenceStart);
                    // Okay, given the alignment will either be at start or at end, we are going to attempt one flip here, might be off slightly
                    int moveBack = deltToUse.FirstSequenceStart != 0L ? 0 - (int)deltToUse.FirstSequenceStart : 0;
                    int start = moveBack + (int)deltToUse.SecondSequenceStart;
					if (start > 0) {
                        var front = _sequence.Skip(start);
                        var back = _sequence.Take(start);
                        var newSeq = new List<byte>(front);
                        newSeq.AddRange(back);
                        var _n = new Sequence(DnaAlphabet.Instance, newSeq.ToArray());
                        if (_n.Count != _sequence.Count) { throw new Exception("Screw up when aligning the sequence"); }
                        if (deltToUse.IsReverseQueryDirection)
                        {
                            ReversePath();
                            _n = _n.GetReverseComplementedSequence() as Sequence;
                            
                        }
                        _sequence = _n;
                    }
                }
                // The underlying bytes should only be accessible from the sequence now
                contigSequence = null;
                finalized = true;
            }
        }




        public int IndexOf(DeBruijnNode item)
        {
            return constituentNodes.IndexOf(item);
        }

        public void Insert(int index, DeBruijnNode item)
        {
            throw new NotImplementedException();
        }

        public void RemoveAt(int index)
        {
            throw new NotImplementedException();
        }

        public DeBruijnNode this[int index]
        {
            get
            {
                return constituentNodes[index];
            }
            set
            {
                throw new NotImplementedException();
            }
        }

        public void Add(DeBruijnNode item)
        {
            throw new NotImplementedException();
        }

        public void Clear()
        {
            throw new NotImplementedException();
        }

        public bool Contains(DeBruijnNode item)
        {
            return nodeHash.Contains(item);
        }

        public void CopyTo(DeBruijnNode[] array, int arrayIndex)
        {
            constituentNodes.CopyTo(array, arrayIndex);
        }

        public int Count
        {
            get { return constituentNodes.Count; }
        }

        public bool IsReadOnly
        {
            get { return true; }
        }

        public bool Remove(DeBruijnNode item)
        {
            throw new NotImplementedException();
        }

        public IEnumerator<DeBruijnNode> GetEnumerator()
        {
            return constituentNodes.GetEnumerator();
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return constituentNodes.GetEnumerator();
        }
    }
}

