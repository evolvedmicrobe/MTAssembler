//#define NO_R
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;
using System.Globalization;
using System.Threading;
using Bio.Algorithms;
using Bio.Algorithms.Assembly;
using Bio.Algorithms.Assembly.Graph;
using Bio.Algorithms.Assembly.Padena;
using Bio.Algorithms.Kmer;
using System.Diagnostics;
using MitoDataAssembler.Visualization;
using System.IO;
using System.Reflection;
using HaploGrepSharp;
using MitoDataAssembler.IndelCaller;
using Bio.Variant;

namespace MitoDataAssembler
{
	//This line of code will take a sequence of reads believed to be similar to the mitochondria
	//it will then create a graph for a de novo assembly
	//Next it will paint that graph with the positions of the kmers in the genome
	//and finally look for candidate junctions where the data doesn't match the genome
	public class MitoPaintedAssembler : ParallelDeNovoAssembler
	{

		#region REPORTVALUES
		/// <summary>
		/// The diagnostic file output prefix, if set a series of diagnostic files are created.
		/// </summary>
		[OutputAttribute]
		public string DiagnosticFileOutputPrefix = null;
		[OutputAttribute]
		public bool SuccessfulAssembly = false;
		[OutputAttribute]
		public string BestMatchingHaplotype = "NOTSET";
		[OutputAttribute]
		public string SecondBestMatchingHaplotype = "NOTSET";
		[OutputAttribute]
		public double BestHaplotypeScore = -1;
		[OutputAttribute]
		public double SecondBestHaplotypeScore = -1;
		[OutputAttribute]
		public int NumberOfEquallyGoodHaplotypes = -1;
		[OutputAttribute]
		public int PolymorphismsMatchingHaplotype = -1;
		[OutputAttribute]
		public int PolymorphismsMissingFromHaplotype = -1;
		[OutputAttribute]
		public int PolymorphismsMissingFromGenotype = -1;
		[OutputAttribute]
		public double PercentNodesPainted;
		[OutputAttribute]
		public double PercentNodesRemovedByLowCoverageOrThreshold;
		//[OutputAttribute]
		//public double MedianKmerCoveragePreFilter;
		//[OutputAttribute]
		//public double MedianKmerCoveragePostFilter;
		[OutputAttribute]
		public long NodeCountAfterCreation;
		[OutputAttribute]
		public long ReferenceNodeCountAfterCreation;
		[OutputAttribute]
		public long NodeCountAfterUndangle;
		[OutputAttribute]
		public long NodeCountAfterRedundancyRemoval;
		[OutputAttribute]
		public long NodeCountAfterCoveragePurge;
		[OutputAttribute]
		public int FinalMetaNodeCount;
		[OutputAttribute]
		public long N50;
		[OutputAttribute]
		public int DecidedAssemblyTotalLength;
		[OutputAttribute]
		public int KmerCutOff;
		[OutputAttribute]
		public long ReadCount;
		[OutputAttribute]
		public long SkippedReadsAfterQCCount;
		[OutputAttribute]
		public int DeletionsFound;
		[OutputAttribute]
		public int PossibleAssemblyCount;
		[OutputAttribute]
		public double PercentageOfScannedReadsUsed;
		[OutputAttribute]
		public long TotalSequencingBP;
		[OutputAttribute]
		public bool DeletionSearchAttempted;
		[OutputAttribute]
		public int SuccessfulAssemblyLength;
		[OutputAttribute]
		public double MinSplitPercentage = -1.0;

        
		#endregion


        /// <summary>
        /// Should we force the sqrt threshold to be used? (Rather than taking the minimum of that or 10?)
        /// </summary>
        public bool ForceSqrtThreshold = false;

        /// <summary>
        /// How far should we go looking for indels?
        /// </summary>
        public int MaxIndelPath = 150;

        /// <summary>
        /// The graph will remove all nodes less than this value for 
        /// before looking for indels if this value is less than the sqrt
        /// of median coverage.
        /// </summary>
        public int AlternateMinimumNodeCount = 0;

        /// <summary>
        /// Output histograms of node counts and graphs before the end of the assembly?
        /// </summary>
        public bool OutputIntermediateGraphSteps = false;

        /// <summary>
        /// Should contig output be skipped?
        /// </summary>
        public bool NoContigOutput = false;

		private HaploGrepSharp.NewSearchMethods.HaploTypeReport haplotypeReport;

		/// <summary>
		/// This is the interface to R, can be used to evaluate R statements.
		/// </summary>
		private RInterface rInt;

		/// <summary>
		/// Should we output diagnostic information?
		/// </summary>
		/// <value><c>true</c> if we output diagnostic information; otherwise, <c>false</c>.</value>
		public bool OutputDiagnosticInformation {
			get{ return !String.IsNullOrWhiteSpace (DiagnosticFileOutputPrefix); }
		}

		/// <summary>
		/// Gets all the nodes in the reference that are derived from 
		/// the rCRS sequence.
		/// </summary>
		/// <value>The reference nodes.</value>
		IEnumerable<DeBruijnNode> referenceNodes {
			get {
				if (Graph != null) {
					return	this.Graph.GetNodes ().Where (x => x.IsInReference);
				} else
					return null;
			}
		}

		public MitoPaintedAssembler () : base ()
		{
			#if !NO_R
			rInt = new RInterface ();
			#endif
		}

		public AssemblyReport GetReport ()
		{
			return new AssemblyReport (this);
		}

		/// <summary>
		/// 
		/// Nigel Rewrite of old method
		/// Estimates and sets erosion and coverage threshold for contigs.
		/// Median value of kmer coverage is set as default value.
		/// Reference: ABySS Release Notes 1.1.1 - "The default threshold 
		/// is the square root of the median k-mer coverage".
		/// </summary>
		protected override void EstimateDefaultThresholds ()
		{

			if (this.AllowErosion || this.AllowLowCoverageContigRemoval) {
                this.EstimateDefaultValuesStarted ();
                // In case of low coverage data, set default as 2.
				// Reference: ABySS Release Notes 1.0.15
				// Before calculating median, discard thresholds less than 2.
				List<long> kmerCoverage = this.Graph.GetNodes ().AsParallel ().Aggregate (
					                                      new List<long> (),
					                                      (kmerList, n) => {
						if (n.KmerCount > 2) {
							kmerList.Add (n.KmerCount);
						}

						return kmerList;
					});

				double threshold;
				if (kmerCoverage.Count == 0) {
					threshold = 2; // For low coverage data, set default as 2
                  
				} else {
					kmerCoverage.Sort ();
					int midPoint = kmerCoverage.Count / 2;
					double median = (kmerCoverage.Count % 2 == 1 || midPoint == 0) ?
                        kmerCoverage [midPoint] :
                        ((float)(kmerCoverage [midPoint] + kmerCoverage [midPoint - 1])) / 2;
					// MedianKmerCoveragePreFilter = median;
					threshold = Math.Sqrt (median);                    
				}
				// Set coverage threshold
				if (this.AllowLowCoverageContigRemoval && this.ContigCoverageThreshold == -1) {
					this.ContigCoverageThreshold = threshold;
				}
				if (this.AllowErosion && this.ErosionThreshold == -1) {
					// Erosion threshold is an int, so round it off
					this.ErosionThreshold = (int)Math.Round (threshold);
				}
                this.EstimateDefaultValuesEnded();				
			}
		}

		private int CalculateSqrtOfMedianCoverageCutoff ()
		{
			var counts = Graph.GetNodes ().Where (x => x.IsInReference).Select (y => y.KmerCount).ToList ();
			counts.Sort ();
			double val = (double)counts [counts.Count / 2];
			val = Math.Sqrt (val);
			return (int)Math.Round (val);
		}

		private void OutputNodeCountHistograms (string FileSuffix, int? cutoff = null)
		{
			List<string> additionalCommands = new List<string> ();
			if (cutoff.HasValue) {
				string cmd = "lines(c(" + cutoff.Value.ToString () + "," + cutoff.Value.ToString () + "),c(0,1000000),lwd=4,col=\"red\")";
				additionalCommands.Add (cmd);
			}
			var refCounts = Graph.GetNodes ().Where (y => y.IsInReference).Select (x => (double)x.KmerCount);
			#if !NO_R
			rInt.HistPDF (refCounts, DiagnosticFileOutputPrefix + @"Refcounts_" + FileSuffix + ".pdf", 125, FileSuffix, "Ref K-mer Occurence", "Count", additionalCommands);
			var allCounts = Graph.GetNodes ().Select (x => (double)x.KmerCount);
			rInt.HistPDF (allCounts, DiagnosticFileOutputPrefix + @"Allcounts_" + FileSuffix + ".pdf", 125, FileSuffix, "All K-mer Occurence", "Count", additionalCommands);
			#endif
		}

		public override PadenaAssembly Assemble (IEnumerable<ISequence> inputSequences)
		{
			if (inputSequences == null) {
				throw new ArgumentNullException ("inputSequences");
			}
			// Step 0: Load the reference genome as a fasta file.
			// Remove ambiguous reads and set up fields for assembler process
			this.Initialize ();

			// Step 1, 2: Create k-mers from reads and build de bruijn graph and paint them with the reference
			System.Diagnostics.Stopwatch sw = System.Diagnostics.Stopwatch.StartNew ();
			this.CreateGraphStarted ();
			this.CreateGraph (inputSequences);
			this.CreateGraphEnded ();
			sw.Stop ();
			this.NodeCountReport ();
			this.TaskTimeSpanReport (sw.Elapsed);           
            
			int count = this.Graph.GetNodes ().Where (x => x.IsInReference).Count ();
			ReferenceNodeCountAfterCreation = count;
			TotalSequencingBP = Graph.GetNodes ().Sum (x => x.KmerCount * KmerLength);
			RaiseMessage ("A total of: " + count.ToString () + " nodes remain from the reference");
			RaiseMessage ("A total of: " + this.Graph.NodeCount + " nodes are in the graph");
            
			NodeCountAfterCreation = Graph.NodeCount;
			SkippedReadsAfterQCCount = Graph.SkippedSequencesCount;
			ReadCount = Graph.ProcessedSequencesCount;

            if (NodeCountAfterCreation < 100) {
				return null;
			}
            
			//Step 2.1, Remove nodes are below coverage cutoff
			sw.Reset ();
			sw.Start ();
            // Estimate and set default value for erosion and coverage thresholds
            this.EstimateDefaultThresholds();
            int sqrtCoverageCutOff = this.CalculateSqrtOfMedianCoverageCutoff();
            var coverageCutOff = ForceSqrtThreshold ? sqrtCoverageCutOff : Math.Min(sqrtCoverageCutOff, AlternateMinimumNodeCount);
            if (MTAssembleArguments.MinNodeCountSet)
            {
                coverageCutOff = AlternateMinimumNodeCount;
            }
            //var coverageCutOff = sqrtCoverageCutOff;
            KmerCutOff = coverageCutOff;
            if (OutputIntermediateGraphSteps)
            {
                OutputNodeCountHistograms("PreFiltered", coverageCutOff);
            }          
			long originalNodes = this.Graph.NodeCount;
			ThresholdCoverageNodePurger snr = new ThresholdCoverageNodePurger (coverageCutOff);
			snr.RemoveLowCoverageNodes (Graph);
			PercentNodesRemovedByLowCoverageOrThreshold = originalNodes / (double)this.Graph.NodeCount;
			sw.Stop ();
            
			TaskTimeSpanReport (sw.Elapsed);
			RaiseMessage ("Finished removing nodes with less than " + snr.CoverageCutOff.ToString () + " counts");
			NodeCountReport ();
			NodeCountAfterCoveragePurge = Graph.NodeCount;
			sw.Reset ();
			sw.Start ();
			RaiseMessage ("Start removing unconnected nodes");

            // Remove pathological nodes
            var badNodes = PathologicalSequencePurger.RemovePathologicalNodes(Graph);
            if (badNodes > 0)
            {
                RaiseMessage("WARNING!!!! Found and removed " + badNodes.ToString()+ " pathological nodes.  These were removed.",false);
            }
                        
            //Step 2.2 Remove nodes that are not connected to the reference genome 
			UnlinkedToReferencePurger remover = new UnlinkedToReferencePurger ();
			remover.RemoveUnconnectedNodes (Graph, referenceNodes);
			RaiseMessage ("Finished removing unconnected nodes");
			this.NodeCountReport ();
			NodeCountAfterUndangle = Graph.NodeCount;
			//outputVisualization ("PostUnconnectedFilter");			           
           
            // Step 2.3: Remove dangling links from graph
			///NIGEL: This also removes the low coverage nodes
			sw.Reset ();
			sw.Restart ();
			this.UndangleGraphStarted ();
			this.UnDangleGraph ();
			this.UndangleGraphEnded ();
			sw.Stop ();
           
            this.TaskTimeSpanReport (sw.Elapsed);
			this.NodeCountReport ();
            if (OutputIntermediateGraphSteps)
            {
                outputVisualization("PostUndangleFilter");
            }

            // Step 4: Remove redundant SNP and indel paths from graph
			RaiseMessage (string.Format (CultureInfo.CurrentCulture, "Starting to remove redundant paths", DateTime.Now));
            this.RemoveRedundancyStarted ();
            RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Starting to remove SNPs", DateTime.Now));
			this.RemoveRedundancy ();
            RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Finished  removing SNPs", DateTime.Now));
            this.NodeCountReport();
            //Now remove redundant indel paths as well
            //TODO: For historic reasons this is largely similar to the snp remover, which isn't so great...
            RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Starting to call INDELs", DateTime.Now));
            var indels = CallAndRemoveIndels();
            
            RaiseMessage(string.Format(CultureInfo.CurrentCulture, "Finished calling and removing small INDELs paths", DateTime.Now));
            this.NodeCountReport();
            // Perform dangling link purger step once more.
            // This is done to remove any links created by redundant paths purger.
            this.UnDangleGraph ();
			RaiseMessage (string.Format (CultureInfo.CurrentCulture, "Finished removing all redundant paths", DateTime.Now));
            this.NodeCountReport();
            
			//STEP 4.2 Rerun the unlinked to reference purger after graph is cleaned
			ChangeNodeVisitFlag (false);
			remover = new UnlinkedToReferencePurger ();
			remover.RemoveUnconnectedNodes (Graph, referenceNodes);
			this.RemoveRedundancyEnded ();
			this.NodeCountReport ();            

			NodeCountAfterRedundancyRemoval = Graph.NodeCount;
            if (OutputIntermediateGraphSteps)
            {
                outputVisualization("Post-redundant-path-removal");
            }
			//Now attempt to assemble and find deletions
			var attemptedAssembly = new MitochondrialAssembly (Graph, DiagnosticFileOutputPrefix);
            FinalMetaNodeCount = attemptedAssembly.AllNodesInGraph.Count;
            SuccessfulAssembly = attemptedAssembly.SuccessfulAssembly;
			if (SuccessfulAssembly) {
				SuccessfulAssemblyLength = attemptedAssembly.AssemblyLength;
				MinSplitPercentage = attemptedAssembly.MinimumGreedySplit;
				if (OutputDiagnosticInformation) {                   
					var outReport = attemptedAssembly.OutputAssembly (DiagnosticFileOutputPrefix);
					if (outReport != null) {
						//TODO: This matching is really crappy, need to find a better way to propogate this on up.
						BestHaplotypeScore = outReport.BestHit.Rank;
						SecondBestHaplotypeScore = outReport.SecondBestHit.Rank;
						BestMatchingHaplotype = outReport.BestHit.node.haplogroup.id;
						SecondBestMatchingHaplotype = outReport.SecondBestHit.node.haplogroup.id;
						NumberOfEquallyGoodHaplotypes = outReport.NumberOfEquallyGoodBestHits;
						PolymorphismsMatchingHaplotype = outReport.BestHit.NumberOfMatchingPolymorphisms;
						PolymorphismsMissingFromHaplotype = outReport.BestHit.NumberOfPolymorphismsMissingFromHaplotype;
						PolymorphismsMissingFromGenotype = outReport.BestHit.NumberOfPolymorphismsMissingFromGenotype;
					}
				} else {
					RaiseMessage ("Greedy assembly skipped as assembly failed.");
				}
			} else {
				RaiseMessage ("Greedy assembly skipped as assembly failed.");               
			}

			//Now find deletions
			this.OutputGraphicAndFindDeletion (attemptedAssembly);            
			PercentageOfScannedReadsUsed = Graph.GetNodes ().Sum (x => x.KmerCount * KmerLength) / (double)TotalSequencingBP;
			RaiseMessage("Used a total of " + PercentageOfScannedReadsUsed.ToString ("p") + " of basepairs in reads for assembly");
            
			// Step 5: Build Contigs - This is essentially independent of deletion finding
			this.BuildContigsStarted ();
			List<ISequence> contigSequences = this.BuildContigs ().ToList ();
			contigSequences.ForEach (x => ReferenceGenome.AssignContigToMTDNA (x));
			contigSequences.Sort ((x, y) => -x.Count.CompareTo (y.Count));
			this.BuildContigsEnded ();

			PadenaAssembly result = new PadenaAssembly ();
			result.AddContigs (contigSequences);
			long totalLength = contigSequences.Sum (x => x.Count);
			RaiseMessage ("Assembled " + totalLength.ToString () + " bases of sequence in " + contigSequences.Count.ToString () + " contigs.");
			if (contigSequences.Count > 0) {                
				N50 = CalculateN50 (contigSequences, totalLength);
				RaiseMessage ("N50: " + N50.ToString ());
			}
			count = this.Graph.GetNodes ().Where (x => x.IsInReference).Count ();
			RaiseMessage ("A total of: " + count.ToString () + " nodes remain from the reference");
			RaiseMessage ("A total of: " + this.Graph.NodeCount + " nodes are in the graph");
			return result;
		}

		protected void OutputGraphicAndFindDeletion (MitochondrialAssembly attemptedAssembly)
		{
			if (attemptedAssembly.SuccessfulAssembly) {
				MitochondrialAssemblyPlotMaker plotMaker = new MitochondrialAssemblyPlotMaker (attemptedAssembly);
				#if !NO_R
                if (OutputIntermediateGraphSteps)
                {
                    plotMaker.Render(rInt, DiagnosticFileOutputPrefix + "AssemblyView.pdf");
                }
				#endif
				DecidedAssemblyTotalLength = plotMaker.Assembly.AssemblyLength;

				//Output all possible assemblies and deletions if possible
				RaiseMessage("Graph contains " + plotMaker.Assembly.AllNodesInGraph.Count.ToString () + " Contained Nodes");
				if (plotMaker.Assembly.AllNodesInGraph.Count < 10) {
					DeletionSearchAttempted = true;
					LargeDeletionFinder ldf = new LargeDeletionFinder ();
					var deletions = ldf.FindAllDeletions (this.Graph, plotMaker.Assembly);
					DeletionsFound = deletions.Count;
					PossibleAssemblyCount = ldf.PossibleDeletionPaths.Count;
					RaiseMessage ("Found a total of: " + deletions.Count + " possible mutations in " + ldf.PossibleDeletionPaths.Count.ToString () + " possible assembly paths");
					//throw error as not finalized here
					ldf.OutputReport (this.DiagnosticFileOutputPrefix + "DeletionReport.csv");
				} else {
					PossibleAssemblyCount = 999;
					DeletionSearchAttempted = false;
				}
			}
		}

		private long CalculateN50 (IList<ISequence> contigs, long total)
		{
			double half = total / 2.0;
			double runningSum = 0;
			foreach (ISequence seq in contigs) {
				runningSum += seq.Count;
				if (runningSum >= half)
					return seq.Count;
			}
			throw new Exception ();
		}

		private int outputVisualization (string prefix)
		{
			if (OutputDiagnosticInformation) {
				OutputNodeCountHistograms (prefix);
				GraphGenerator gg = new GraphGenerator (this.Graph);
				gg.OutputGraph (DiagnosticFileOutputPrefix + prefix + ".graphml");
				if (gg.MetaNodes != null) {
					return gg.MetaNodes.Count;
				}
			}
			return 0; 
		}

		/// <summary>
		/// Step 1: Building k-mers from sequence reads
		/// Step 2: Build de bruijn graph for input set of k-mers.
		/// Sets the _assemblerGraph field.
		/// </summary>
        protected override void CreateGraph(IEnumerable<ISequence> inputSequences )
		{
			this.Graph = new DeBruijnGraph (this._kmerLength);
			this.Graph.Build (inputSequences, false);
			RaiseMessage ("Graph Processed:\t" + this.Graph.ProcessedSequencesCount.ToString () + " sequences");
			RaiseMessage ("Skipped:\t" + this.Graph.SkippedSequencesCount.ToString () + " sequences");
			// Recapture the kmer length to keep them in sync.
			this._kmerLength = this.Graph.KmerLength;
			PaintKmersWithReference ();           
			this.Graph.DestroyKmerManager ();
		}

		private void OutCSV (Dictionary<long, long> histogram, string name)
		{
			System.IO.StreamWriter SW = new System.IO.StreamWriter (name);
			SW.WriteLine ("KmerCoverage,Frequency");
			foreach (var v in histogram.Select(x => new { coverage = x.Key, freq = x.Value }).OrderBy(y => y.coverage)) {
				SW.WriteLine (v.coverage.ToString () + "," + v.freq.ToString ());
			}
			SW.Close ();
		}

		/// <summary>
		/// Changes the visit flag on all nodes so it is false.
		/// </summary>
		/// <param name="valueToSet">If set to <c>true</c> value to set.</param>
		public void ChangeNodeVisitFlag (bool valueToSet)
		{
			Parallel.ForEach (this.Graph.GetNodes (), node => node.IsVisited = valueToSet);
		}

		/// <summary>
		/// Add a line to each debruijin node if it corresponds to a 
		/// kmer from a single position in a reference genome, 
		/// </summary>
		protected void PaintKmersWithReference ()
		{
			List<int> missingLocs = new List<int> ();
			var refKmerPositions = SequenceToKmerBuilder.BuildKmerDictionary (ReferenceGenome.ReferenceSequence, this.KmerLength);
			int KmersPainted = 0;
			int KmersSkipped = 0;
			DeBruijnGraph graph = this.Graph;
			long totalNodes = graph.NodeCount;
			foreach (var v in refKmerPositions) {
				ISequence seq = v.Key;
				IList<long> locations = v.Value;
				if (locations.Count == 1) {
					var kmerData = new KmerData32 ();
					kmerData.SetKmerData (seq, 0, this.KmerLength);
					DeBruijnNode matchingNode = this.Graph.KmerManager.SetNewOrGetOld (kmerData, false);
					if (matchingNode != null) {
						matchingNode.ReferenceGenomePosition = (short)locations [0];
						KmersPainted++;
						if (matchingNode.ReferenceGenomePosition < 0)
							throw new Exception ();
					} else {
						missingLocs.Add ((int)locations [0]);
					}
				} else {
					KmersSkipped += locations.Count;
				}
			}
            if (false && OutputDiagnosticInformation)
            {
                StreamWriter sw = new StreamWriter("OutMissing.csv");
                foreach (int i in missingLocs)
                {
                    sw.WriteLine(i.ToString());
                }
                sw.Close();
            }
			double percentKmersSkipped = 100.0 * (KmersSkipped) / ((double)(KmersPainted + KmersSkipped));
			if (percentKmersSkipped > 95.0) {
				throw new InvalidProgramException ("Reference Genome Skipped over 95% of Kmers");
			}
			double percentHit = KmersPainted / (double)refKmerPositions.Count;
			RaiseMessage ("A total of " + (100.0 * percentHit).ToString () + "% nodes in the reference were painted");
			PercentNodesPainted = 100.0 * KmersPainted / (double)totalNodes;
			RaiseMessage(PercentNodesPainted.ToString ("n2") + " % of nodes painted, for a total of " + KmersPainted.ToString () + " painted.");
            RaiseMessage(percentKmersSkipped.ToString ("n2") + " % of Kmers were skipped for being in multiple locations");
		}

		/// <summary>
		/// Step 5: Build contigs from de bruijn graph.
		/// If coverage threshold is set, remove low coverage contigs.
        /// 
		/// </summary>
		/// <returns>List of contig sequences.</returns>
		protected override IEnumerable<ISequence> BuildContigs ()
		{
			if (this.ContigBuilder == null) {
				throw new Exception ();
			}
			// Step 5.1: Remove low coverage contigs
			if (this.AllowLowCoverageContigRemoval && this.ContigCoverageThreshold > 0) {
				this.LowCoverageContigPurger.RemoveLowCoverageContigs (this.Graph, this.ContigCoverageThreshold);
			}

			// Step 5.2: Build Contigs
			return this.ContigBuilder.Build (this.Graph);
		}

        /// <summary>
        /// Calls indels up to a certain depth and removes redundant paths except those with the highest k-mer coverage.
        /// </summary>
        /// <returns></returns>
        protected List<ContinuousFrequencyIndelGenotype> CallAndRemoveIndels()
        {
            //DeBruijnPathList redundantNodes;
            try
            {
                ContinuousGenotypeIndelCaller indelCaller = new ContinuousGenotypeIndelCaller(MaxIndelPath);
                return indelCaller.CallAndRemoveIndels(this.Graph);
            }
            catch (Exception thrown)
            {
                throw new Exception(
 @"Program failed attempting to call Indels.
A likely cause is that there were so many nodes in the graph that it ran out of memory. Please check that the total nodes are near ~16,569 by this step, if the node count is high re-running the program with the -force_sqrt option may solve this.

The exact error returned was:
" + thrown.Message, thrown);  
            }
        }
	}
}