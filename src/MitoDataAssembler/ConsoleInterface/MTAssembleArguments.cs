using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using Bio;
using Bio.Algorithms.Assembly;
using Bio.Algorithms.Assembly.Padena;
using Bio.IO.BAM;
using Bio.Util;
using Bio.IO.FastA;
using MitoDataAssembler.PairedEnd;

namespace MitoDataAssembler
{
	/// <summary>
	/// Class for Assemble options.  Kind of a misnomer because it also does the assembly...
    /// TODO: Factor out options from actions on them
	/// </summary>
	public class MTAssembleArguments
	{

		#region Public Fields
        /// <summary>
		/// Length of k-mer.
		/// </summary>
		public int KmerLength = 19;
		/// <summary>
		/// Threshold for removing dangling ends in graph.
		/// </summary>
		public int DangleThreshold = -1;
		/// <summary>
		/// Length Threshold for removing redundant paths in graph.
		/// </summary>
		public int RedundantPathLengthThreshold = -1;
		/// <summary>
		/// Threshold for eroding low coverage ends.
		/// </summary>
		public int ErosionThreshold = -1;
		/// <summary>
		/// Bool to do erosion or not.
		/// </summary>
		public bool AllowErosion = false;


        /// <summary>
        /// The minimum number of times a k-mer must appear before it is considered for the 
        /// assembly and indel finding.  (or the Sqrt of the median coverage).
        /// </summary>
        public int MinimumNodeCount = 11;

        /// <summary>
        /// HACK: This is a value to check if the user has set the node count, in which case WE MUST go with it!
        /// </summary>
        public static bool MinNodeCountSet = false;

        /// <summary>
        /// Should we force the sqrt threshold to be used? (Rather than taking the minimum of that or 10?)
        /// </summary>
        public bool ForceSqrtThreshold = false;

        /// <summary>
        /// Should contig output be skipped?
        /// </summary>
        public bool NoContigOutput = false;

		/// <summary>
		/// Whether to estimate kmer length.
		/// </summary>
		public bool AllowKmerLengthEstimation = false;
		/// <summary>
		/// Threshold used for removing low-coverage contigs.
		/// </summary>
		public int ContigCoverageThreshold = -1;

		/// <summary>
		/// Help.
		/// </summary>
		public bool Help = false;
		/// <summary>
		/// Input file of reads.
		/// </summary>
		public string Filename = string.Empty;
		private string fullFileName;

        /// <summary>
		/// Skip calling SNPs and Haplotypes using a column wise pile-up?
        /// </summary>
		public bool Skip_Pileup_Calling = false;

        /// <summary>
        /// Are we going to use the EM algorithm to do frequency estimates, 
        /// or just use read counts from the pile-ups??
        /// </summary>
        public bool Skip_EM_Frequency_Estimates = false;

		/// <summary>
		/// Should we skip the denovo assembler?
		/// </summary>
		public bool Skip_Assembly_Step = false;


		/// <summary>
		/// Should we skip the peak finding step?
		/// </summary>
		public bool Skip_Peak_Finder = false;

		/// <summary>
		/// If a subset of the BAM is used this can be specified here.
		/// </summary>
		public string ChromosomeName = "MT";

		/// <summary>
		/// Prefix of the report output prefix file
		/// </summary>
		public string ReportOutputSuffix = "Report";

        /// <summary>
        /// Make a depth of coverage plot.
        /// </summary>
        public bool Skip_DepthOfCoveragePlot = false;


        /// <summary>
        /// Output histograms of node counts and graphs before the end of the assembly?
        /// </summary>
        public bool OutputIntermediateGraphSteps = false;

		/// <summary>
		/// Display verbose logging during processing.
		/// </summary>
		public bool Verbose {
			get { return _pVerbose; }
			set {
				if (value)
					Output.TraceLevel = OutputLevel.Information | OutputLevel.Verbose;
				else
					Output.TraceLevel = OutputLevel.Information;
				_pVerbose = value;
				//TODO: Make this an option again
				_pVerbose = true;

			}
		}

		private bool _pVerbose = true;
		/// <summary>
		/// Quiet flag (no logging)
		/// </summary>
		public bool Quiet = false;
		/// <summary>
		/// The diagnostic file prefix, will appear in front of all output.
		/// </summary>
		public string DiagnosticFilePrefix = string.Empty;


		public string ContigFileName { get { return DiagnosticFilePrefix + "contigs.fna"; } }
		#endregion

		#region Public methods
        public void OutputEnvironmentSettings()
        {
            Output.WriteLine(OutputLevel.Verbose, "\nEnvironment Settings:");
            Output.WriteLine(OutputLevel.Verbose, "\tk-mer Length: {0}", this.KmerLength);
            Output.WriteLine(OutputLevel.Verbose, "\tPrefix is: " + DiagnosticFilePrefix);
            Output.WriteLine(OutputLevel.Verbose, "\tDiagnostic Information On: " + (!Quiet).ToString());


            string dll = System.Environment.GetEnvironmentVariable( RInterface.R_LIB_ENV_DIR);
            string r_home = System.Environment.GetEnvironmentVariable(RInterface.R_HOME_ENV_DIR);
            Output.WriteLine(OutputLevel.Verbose, "\t"+RInterface.R_HOME_ENV_DIR + ": " + r_home);
            Output.WriteLine(OutputLevel.Verbose, "\t"+RInterface.R_LIB_ENV_DIR + ": " + dll);
            var platform  = RDotNet.NativeLibrary.NativeUtility.GetPlatform();
            Output.WriteLine(OutputLevel.Verbose, "\tPlatform: "+platform.ToString());

        }

		/// <summary>
		/// Assembles the sequences and returns the string that can be placed in a CSV output report.
		/// </summary>
		/// <returns></returns>
		public List<AlgorithmReport> ProcessMTDNA()
		{
			var results = new List<AlgorithmReport> ();

			//STEP 0: Preprocess reads file 
			FileInfo refFileinfo = new FileInfo (this.Filename);
			long refFileLength = refFileinfo.Length;


			fullFileName = Path.GetFullPath (this.Filename);
			if (File.Exists (fullFileName)) {
				Output.WriteLine (OutputLevel.Verbose);
				Output.WriteLine (OutputLevel.Verbose, "Found read file: {0}", fullFileName);
				Output.WriteLine (OutputLevel.Verbose, "   File Size           : {0}", Program.FormatMemorySize(refFileLength));
			}
			//Print environment settings to console.  
            OutputEnvironmentSettings();
            
			//Assemble
			var asmReport = CreateAssemblyAndDepthOfCoverage ();
			results.Add (asmReport);

			//Peak find
			var peakFindReport = RunPeakFinder ();
			results.Add (peakFindReport);

            //Pile-up
			var pileupReport = RunSNPCaller ();
			results.Add (pileupReport);

			if (!String.IsNullOrEmpty (DiagnosticFilePrefix)) {
                var outFile = new StreamWriter(DiagnosticFilePrefix + ReportOutputSuffix + ".csv");
				var header = String.Join (",", results.Select (z => z.HeaderLineForCSV ));
				outFile.WriteLine (header);
				var data = String.Join (",", results.Select (z => z.DataLineForCSV ));
				outFile.WriteLine (data);
				outFile.Close ();
			}
			return results;
		}

		#endregion

		#region Protected Members

		protected AssemblyReport CreateAssemblyAndDepthOfCoverage() 
		{
			if (Skip_Assembly_Step) {
				return new AssemblyReport ();
			}

            DepthOfCoverageGraphMaker coveragePlotter = !Skip_DepthOfCoveragePlot ?
                                                        new DepthOfCoverageGraphMaker() : null;

            IEnumerable<ISequence> reads = this.createSequenceProducer(this.Filename, coveragePlotter , true);
            TimeSpan algorithmSpan = new TimeSpan();
            Stopwatch runAlgorithm = new Stopwatch();
		

            //Step 1: Initialize assembler.	
            Output.WriteLine(OutputLevel.Verbose, "\nAssemblying mtDNA and obtaining depth of coverage (if asked).");
            MitoPaintedAssembler.StatusChanged += this.StatusChanged;
            MitoPaintedAssembler assembler = new MitoPaintedAssembler() { 
                DiagnosticFileOutputPrefix = DiagnosticFilePrefix,
                AllowErosion = AllowErosion,
                AlternateMinimumNodeCount = MinimumNodeCount,
                DanglingLinksThreshold = DangleThreshold,
                ErosionThreshold = ErosionThreshold,
                AllowKmerLengthEstimation = AllowKmerLengthEstimation,
                RedundantPathLengthThreshold = RedundantPathLengthThreshold,
                OutputIntermediateGraphSteps = OutputIntermediateGraphSteps,
                NoContigOutput = NoContigOutput,
                ForceSqrtThreshold = ForceSqrtThreshold
            };
            if (ContigCoverageThreshold != -1)
            {
                assembler.AllowLowCoverageContigRemoval = true;
                assembler.ContigCoverageThreshold = ContigCoverageThreshold;
            }
            if (!this.AllowKmerLengthEstimation)
            {
                assembler.KmerLength = this.KmerLength;
            }
            
			//Step 2: Assemble
			runAlgorithm.Restart ();
			var assembly = assembler.Assemble (reads);
			runAlgorithm.Stop ();
			algorithmSpan = algorithmSpan.Add (runAlgorithm.Elapsed);
			if (this.Verbose) {
				Output.WriteLine (OutputLevel.Verbose);
				Output.WriteLine (OutputLevel.Verbose, "\tCompute time: {0}", runAlgorithm.Elapsed);
			}

			//Step 3: Report
            if (!NoContigOutput)
            {
                runAlgorithm.Restart();
                this.writeContigs(assembly);
                runAlgorithm.Stop();
            }
			algorithmSpan = algorithmSpan.Add (runAlgorithm.Elapsed);

			if (this.Verbose) {
				Output.WriteLine (OutputLevel.Verbose);
				Output.WriteLine (OutputLevel.Verbose, "\tWrite contigs time: {0}", runAlgorithm.Elapsed);
				Output.WriteLine (OutputLevel.Verbose, "\tTotal assembly runtime: {0}", algorithmSpan);
			}

			if (coveragePlotter !=null) {
				coveragePlotter.OutputCoverageGraphAndCSV (DiagnosticFilePrefix);
			}

			return assembler.GetReport ();
		}

		protected PairedEndDeletionFinderReport RunPeakFinder()
		{
			if (Skip_Peak_Finder) {
				return new PairedEndDeletionFinderReport ();
			}
			//Step 4: Run Break Finder
			Output.WriteLine (OutputLevel.Verbose, "\nRunning Peak Finder Program");
            Stopwatch time = new Stopwatch();
            time.Start();
			var isBAM = Helper.IsBAM (fullFileName);
			PairedEndDeletionFinderReport report;
			try {
				if (isBAM) {
					PairedEndPeakFinder pairedEndPeakFinder = new PairedEndPeakFinder (fullFileName, DiagnosticFilePrefix);
					pairedEndPeakFinder.FindDeletionPeaks ();
					report = new PairedEndDeletionFinderReport (pairedEndPeakFinder);
				} else {
					report = new PairedEndDeletionFinderReport ("NOT BAM");
				}
			} catch (Exception thrown) {
				Output.WriteLine (OutputLevel.Error, "\tFailed to run peak finder: " + thrown.Message);
				report = new PairedEndDeletionFinderReport ("Failure");
			}
            time.Stop();
            if (this.Verbose)
            {
                Output.WriteLine(OutputLevel.Verbose, "\tTime Elapsed: {0}", time.Elapsed);
            }
			return report;
		}

		protected SNPCallerReport RunSNPCaller()
		{
			if (Skip_Pileup_Calling) {
				return new SNPCallerReport (AlgorithmResult.NotAttempted);
			}
            Output.WriteLine(OutputLevel.Verbose, "\nCalling Pile-up SNPs");
			Stopwatch time = new Stopwatch ();
			time.Start ();
			var reads = this.createSequenceProducer (this.Filename);
			Bio.Variant.ContinuousGenotypeCaller.DO_EM_ESTIMATION = !Skip_EM_Frequency_Estimates;
			var res = SNPCaller.CallSNPs (reads );
			time.Stop ();
			if (this.Verbose) {
				Output.WriteLine (OutputLevel.Verbose, "\tTime Elapsed: {0}", time.Elapsed);
			}
			return res;
		}

		/// <summary>
		/// Writes the contigs to the file.
		/// </summary>
		/// <param name="assembly">IDeNovoAssembly parameter is the result of running De Novo Assembly on a set of two or more sequences. </param>
		protected void writeContigs (PadenaAssembly assembly)
		{
			if (assembly.AssembledSequences.Count == 0) {
				Output.WriteLine (OutputLevel.Results, "\tNo sequences assembled.");
				return;
			}
			ensureContigNames (assembly.AssembledSequences);

			if (!string.IsNullOrEmpty (this.DiagnosticFilePrefix)) {
				using (FastAFormatter formatter = new FastAFormatter (ContigFileName)) {
					formatter.AutoFlush = true;
					foreach (ISequence seq in assembly.AssembledSequences) {
						formatter.Write (seq);
					}
				}
				Output.WriteLine (OutputLevel.Information, "\tWrote {0} sequences to {1}", assembly.AssembledSequences.Count, ContigFileName);
			} else {
				Output.WriteLine (OutputLevel.Information, "\tAssembled Sequence Results: {0} sequences", assembly.AssembledSequences.Count);
				using (FastAFormatter formatter = new FastAFormatter ()) {
					formatter.Open (new StreamWriter (Console.OpenStandardOutput ()));
					formatter.MaxSymbolsAllowedPerLine = decideOutputWidth ();
					formatter.AutoFlush = true;
					foreach (ISequence seq in assembly.AssembledSequences) {
						formatter.Write (seq);
					}
				}
			}
		}

		private int decideOutputWidth ()
		{
			if (Console.WindowWidth < 50)
				return 80;
			else
				return Math.Min (80, Console.WindowWidth - 2);
		}

		/// <summary>
		/// Ensures the sequence contigs have a valid ID. If no ID is present
		/// then one is generated from the index and filename.
		/// </summary>
		/// <param name="sequences"></param>
		private void ensureContigNames (IList<ISequence> sequences)
		{
			for (int index = 0; index < sequences.Count; index++) {
				ISequence inputSequence = sequences [index];
				if (string.IsNullOrEmpty (inputSequence.ID))
					inputSequence.ID = generateSequenceId (index + 1);
				else
					inputSequence.ID = generateSequenceId (index + 1) + inputSequence.ID;
			}
		}

		/// <summary>
		/// Generates a sequence Id using the output filename, or first input file.
		/// </summary>
		/// <param name="counter">Sequence counter</param>
		/// <returns>Auto-generated sequence id</returns>
		private string generateSequenceId (int counter)
		{
			string filename = Path.GetFileNameWithoutExtension (this.ContigFileName);
			if (string.IsNullOrEmpty (filename))
				filename = Path.GetFileNameWithoutExtension (this.Filename);
			filename = filename.Replace (" ", "");
			Debug.Assert (!string.IsNullOrEmpty (filename));
			return string.Format (CultureInfo.InvariantCulture, "{0}_{1}", filename, counter);
		}

		
		/// <summary>
		/// Create a sequence enumerator that filters the reads and adds them to the depth of coverage counter
        /// if necessary.
		/// </summary>
		/// <param name="fileName">Filename to load data from</param>
		/// <returns>Enumerable set of ISequence elements</returns>
		private IEnumerable<CompactSAMSequence> createSequenceProducer (string fileName, DepthOfCoverageGraphMaker coveragePlotter = null, bool alsoGetNuclearHits = false)
		{
            if (!Skip_DepthOfCoveragePlot && !Helper.IsBAM(fileName))
            {
                Skip_DepthOfCoveragePlot = true;
                Output.WriteLine(OutputLevel.Error, "Warning: No coverage plots can be made without an input BAM File");
            }
            IEnumerable<CompactSAMSequence> sequences;
            if (!alsoGetNuclearHits)
            {
                var parser = new BAMSequenceParser(Filename);
                if (ChromosomeName != string.Empty)
                {
                    parser.ChromosomeToGet = ChromosomeName;
                    sequences = parser.Parse();
                }
                else
                {
                    sequences = parser.Parse();
                }
            }
            else
            {
                sequences = BAMNuclearChromosomeReadGenerator.GetNuclearAndMitochondrialReads(fileName);
            }
            //Filter by quality
            return ReadFilter.FilterReads(sequences,coveragePlotter);
		}

		/// <summary>
		/// Method to handle status changed event.
		/// </summary>
		protected void StatusChanged (object sender, StatusChangedEventArgs statusEventArgs)
		{
			if (Verbose)
				Output.WriteLine (OutputLevel.Verbose, statusEventArgs.StatusMessage);
			else if (!Quiet) {
				if (statusEventArgs.StatusMessage.StartsWith ("Step", StringComparison.OrdinalIgnoreCase)
				    && statusEventArgs.StatusMessage.Contains ("Start")) {
					int pos = statusEventArgs.StatusMessage.IndexOf (" - ", StringComparison.OrdinalIgnoreCase);
					Output.WriteLine (OutputLevel.Information, statusEventArgs.StatusMessage.Substring (0, pos));
				}
			}
		}



		#endregion
	}
	
}
