using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly.Padena;
using Bio.Algorithms.Assembly.Graph;
using MitoDataAssembler.Visualization;
using MitoDataAssembler.Extensions;
using Bio.IO.FastA;
using System.IO;

namespace MitoDataAssembler
{
    /// <summary>
    /// A mitochondrial assembly, formed from a list of 
    /// </summary>
    public class MitochondrialAssembly
    {
        /// <summary>
        /// A class to hold information about possible splits at a junction that were/could have been followed on the way to the assembly.
        /// </summary>
        public class SplitData
        {
            //The frequency of the highest split reads relative to total reads at this position.
            public double MaxFrequency;
            /// <summary>
            /// A list of possibilities at this junction.
            /// </summary>
            public List<MetaNode.PossiblePath> PossiblePaths;
            /// <summary>
            /// The K-mer count of the most frequent route
            /// </summary>
            public int MaxCount;
            /// <summary>
            /// How many options were there?
            /// </summary>
            public int TotalPathsPossible { get { return PossiblePaths.Count; } }
            public SplitData(List<MetaNode.PossiblePath> paths)
            {
                this.PossiblePaths = paths;
                paths.Sort((x, y) => -x.Weight.CompareTo(y.Weight));
                paths.ForEach(x => { if (x.NeedsReversing) x.NeighborNode.ReversePath(); x.NeedsReversing = false; });
                double sum = paths.Sum(x => x.Weight);
                MaxFrequency = paths[0].Weight / (double)sum;
            }
            public MetaNode.PossiblePath BestPath
            {
                get { return PossiblePaths[0]; }
            }

        }
        
		/// <summary>
		/// The lowest frequency selected along a greedy path
		/// </summary>
		public double MinimumGreedySplit=Double.MaxValue;

        private PossibleAssembly _greedyPathAssembly;
        public PossibleAssembly GreedyPathAssembly
        {
            get
            {
                if (SuccessfulAssembly)
                    return _greedyPathAssembly;
                else
                    return null;
            }
        }
		/// <summary>
		/// The sequence of the forward primer used by the LR-PCR at the Broad and at Baylor
		/// </summary>
		public string forwardPrimer="CCGCACAAGAGTGCTACTCTCCTC";
		/// <summary>
		/// The reverse complement sequence of the reverse primer used by the LR-PCR at the Broad and at Baylor
		/// </summary>
		public string reversePrimer="CACCATCCTCCGTGAAATCAATATC";//Actually the RC of the primer
		public int AssemblyLength=StaticResources.CRS_LENGTH;
        public List<SplitData> PathSplits = new List<SplitData>();
        private GraphGenerator gg;
        public bool FormsCompleteLoop;
        public bool SuccessfulAssembly { get; set; }
        List<MetaNode> assemblyNodes=new List<MetaNode>();
        public List<MetaNode> NodesInCompleteAssembly
        {
            get
            {
                if (SuccessfulAssembly) { return assemblyNodes; }
                else return null;
            }
        }
        public List<MetaNode> AllNodesInGraph
        {
            get { return gg.MetaNodes; }
        }
        public MitochondrialAssembly(DeBruijnGraph graph,string Prefix)
        {
            gg = new GraphGenerator(graph);
            gg.OutputDotGraph(Prefix+"_Graph.dot");
            if (gg.MetaNodes.Count > 0)
            {
                attemptToCreateAssembly();
            }
            else { SuccessfulAssembly = false; FormsCompleteLoop = false; }
        }
        public HaploGrepSharp.NewSearchMethods.HaploTypeReport OutputAssembly(string fileNamePrefix)
        {
			if (SuccessfulAssembly) {
				FastAFormatter fa = new FastAFormatter (fileNamePrefix + "_BestGreedyAssembly.fna");
				StringBuilder sb = new StringBuilder (StaticResources.CRS_LENGTH);
				var bestAssembly = GreedyPathAssembly;
				bestAssembly.Finalize ();
				Bio.Sequence s = new Bio.Sequence (bestAssembly.Sequence);
				s.ID = "GreedyAssembly - length=" + AssemblyLength.ToString ();// + bestAssembly.FirstReferencePosition.Value.ToString() + " - " + GreedyPathAssembly.LastReferencePosition.Value.ToString();
				fa.Write (s);
				fa.Close ();
				//Now report all differences as well
				StreamWriter sw = new StreamWriter (fileNamePrefix + "_Report.txt");
				var searcher = new HaploGrepSharp.NewSearchMethods.HaplotypeSearcher ();
				List<string> linesToWrite = new List<string> ();
				var report = searcher.GetHaplotypeReport (s, linesToWrite, fileNamePrefix);
				foreach (var l in linesToWrite) {
					sw.WriteLine (l);
				}
				sw.Close ();
				return report;
			}
            return null;
        }
        private void attemptToCreateAssembly()
        {
            //TODO: This node should always be a good start node, but may be an erroneous one, check for this.
            var curNode = gg.MetaNodes.Where(x => x.Lowest_Reference_Position != 0).MaxBy(x => (x.AvgKmerCoverage * x.ConstituentNodes.Count));//*(.2/x.Lowest_Reference_Position));//.MinBy(x => x.Lowest_Reference_Position);

            //Let's try just going with the forward primer
            //var match = forwardPrimer.Substring(0, gg.MegaNodes.First().LeadingKmer.Length);
            //var rc_match = ((new Bio.Sequence(Bio.Alphabets.NoGapDNA, match)).GetReverseComplementedSequence() as Bio.Sequence).ConvertToString();
            //var curNode = gg.MegaNodes.Where(x => x.Sequence.Contains(match) || x.Sequence.Contains(rc_match)).First();
            _greedyPathAssembly = new PossibleAssembly();
            if (!curNode.CircularLoop)
            {
                Console.WriteLine("Attempting to find greedy path, frequencies of majority split below");
                //now to attempt to loop back to the start node
                //will move along while greedily grabbing the next node with the highest kmer coverage
                //constantly oriented everyone so we go right ot left
                while (true)
                {
                    assemblyNodes.Add(curNode);
                    _greedyPathAssembly.AddMetaNode(curNode);
                    var possibles = curNode.GetOutgoingNodes().ToList();
                    if (possibles.Count > 0)
                    {
                        SplitData sd = new SplitData(possibles);
                        PathSplits.Add(sd);
                        Console.WriteLine("Possible Paths: "+possibles.Count+"\tFrequency: "+sd.MaxFrequency.ToString());
                        if (sd.MaxFrequency < MinimumGreedySplit)
                        {
                            MinimumGreedySplit = sd.MaxFrequency;
                        }
                        curNode = sd.BestPath.NeighborNode;
                        if (assemblyNodes.Contains(curNode))
                        {
                            FormsCompleteLoop = true;
                            break;
                        }

                    }
                    else { FormsCompleteLoop = false; SuccessfulAssembly = false; break; }
                }
            }
            else
            {
                FormsCompleteLoop = true;
                assemblyNodes.Add(curNode);
                _greedyPathAssembly.AddMetaNode(curNode);
                MinimumGreedySplit = 1.0;
            }
            int length = assemblyNodes.Sum(x => x.LengthOfNode);
            //now, did we form an assembly?
            if (FormsCompleteLoop || Math.Abs(length - AssemblyLength) < 100)
            {
                SuccessfulAssembly = true;
                _greedyPathAssembly.Finalize();
                AssemblyLength = (int) _greedyPathAssembly.Sequence.Count;
                //TODO: More sophisticated criteria than larger than 16 kb to validate assembly
                if (AssemblyLength > 8000)
                {
                    SuccessfulAssembly = true;
                    Console.WriteLine("Successful assembly of length: " + AssemblyLength.ToString());
                }
                else
                {
                    SuccessfulAssembly = false;
                    Console.WriteLine("Assembly failed.  Only recovered sequence of length: " + AssemblyLength.ToString());                   
                }
            }
        }
    }
}
