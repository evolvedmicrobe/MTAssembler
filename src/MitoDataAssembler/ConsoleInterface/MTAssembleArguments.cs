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
		/// Whether to estimate kmer length.
		/// </summary>
		public bool AllowKmerLengthEstimation = false;
		/// <summary>
		/// Threshold used for removing low-coverage contigs.
		/// </summary>
		public int ContigCoverageThreshold = -1;
		/// <summary>
		/// Force specified kmer (no warning prompt)
		/// </summary>
		public bool ForceKmer = false;
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
        /// Should we also call SNPs and Haplotypes using a column wise pile-up?
        /// </summary>
        public bool DoPileUpSNPCalling = true;
		/// <summary>
		/// Output file.
		/// </summary>
		public string OutputFile = string.Empty;

		/// <summary>
		/// If a subset of the BAM is used this can be specified here.
		/// </summary>
		public string ChromosomeName = "MT";

		/// <summary>
		/// Prefix of the report output prefix file
		/// </summary>
		public string ReportOutputPrefix = "AssemblyReport_";

        /// <summary>
        /// Make a depth of coverage plot.
        /// </summary>
        public bool MakeDepthOfCoveragePlot = true;

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

		#endregion

		#region Public methods
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

			Output.WriteLine (OutputLevel.Information, StaticResources.AssemblyStarting);
			fullFileName = Path.GetFullPath (this.Filename);
			if (File.Exists (fullFileName)) {
				Output.WriteLine (OutputLevel.Verbose);
				Output.WriteLine (OutputLevel.Verbose, "Found read file: {0}", fullFileName);
				Output.WriteLine (OutputLevel.Verbose, "   File Size           : {0}", refFileLength);
				Output.WriteLine (OutputLevel.Verbose, "   k-mer Length        : {0}", this.KmerLength);
			}

			Console.WriteLine ("Prefix is: " + DiagnosticFilePrefix);
			Console.WriteLine ("Diagnostic Information On: " + (!Quiet).ToString());

			//Assemble
			var asmReport = CreateAssemblyAndDepthOfCoverage ();
			results.Add (asmReport);

			//Peak find
			var peakFindReport = RunPeakFinder ();
			results.Add (peakFindReport);

            //Pile-up
			if (DoPileUpSNPCalling) {
				var pileupReport = RunSNPCaller ();
				results.Add (pileupReport);
			}


			if (!String.IsNullOrEmpty (DiagnosticFilePrefix)) {
				var outFile = new StreamWriter (ReportOutputPrefix + DiagnosticFilePrefix + ".csv");
				var header = String.Join (",", results.Select (z => z.HeaderLineForCSV ));
				outFile.Write (header);
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

			MitoPaintedAssembler assembler = new MitoPaintedAssembler ();
			DepthOfCoverageGraphMaker coveragePlotter = MakeDepthOfCoveragePlot ? 
														new DepthOfCoverageGraphMaker() : null;
		
			IEnumerable<ISequence> reads = this.createSequenceProducer (this.Filename,coveragePlotter);

			TimeSpan algorithmSpan = new TimeSpan ();
			Stopwatch runAlgorithm = new Stopwatch ();

			//Step 1: Initialize assembler.
			assembler.DiagnosticFileOutputPrefix = DiagnosticFilePrefix;				
			assembler.StatusChanged += this.StatusChanged;
			assembler.AllowErosion = this.AllowErosion;
			assembler.AllowKmerLengthEstimation = this.AllowKmerLengthEstimation;
			if (ContigCoverageThreshold != -1) {
				assembler.AllowLowCoverageContigRemoval = true;
				assembler.ContigCoverageThreshold = ContigCoverageThreshold;
			}
			assembler.DanglingLinksThreshold = this.DangleThreshold;
			assembler.ErosionThreshold = this.ErosionThreshold;
			if (!this.AllowKmerLengthEstimation) {
				assembler.KmerLength = this.KmerLength;
			}
			assembler.RedundantPathLengthThreshold = this.RedundantPathLengthThreshold;

			//Step 2: Assemble
			runAlgorithm.Restart ();
			IDeNovoAssembly assembly = assembler.Assemble (reads);
			runAlgorithm.Stop ();
			algorithmSpan = algorithmSpan.Add (runAlgorithm.Elapsed);
			if (this.Verbose) {
				Output.WriteLine (OutputLevel.Verbose);
				Output.WriteLine (OutputLevel.Verbose, "Compute time: {0}", runAlgorithm.Elapsed);
			}

			//Step 3: Report
			runAlgorithm.Restart ();
			this.writeContigs (assembly);
			runAlgorithm.Stop ();
			algorithmSpan = algorithmSpan.Add (runAlgorithm.Elapsed);

			if (this.Verbose) {
				Output.WriteLine (OutputLevel.Verbose);
				Output.WriteLine (OutputLevel.Verbose, "Write contigs time: {0}", runAlgorithm.Elapsed);
				Output.WriteLine (OutputLevel.Verbose, "Total assembly runtime: {0}", algorithmSpan);
			}

			if (coveragePlotter !=null) {
				coveragePlotter.OutputCoverageGraphAndCSV (DiagnosticFilePrefix);
			}

			return assembler.GetReport ();
		}

		protected PairedEndDeletionFinderReport RunPeakFinder()
		{
			//Step 4: Run Break Finder
			Output.WriteLine (OutputLevel.Verbose, "Attempting to Run Peak Finder Program");
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
				Output.WriteLine (OutputLevel.Error, "Failed to run peak finder: " + thrown.Message);
				report = new PairedEndDeletionFinderReport ("Failure");
			}
			return report;
		}

		protected SNPCallerReport RunSNPCaller()
		{

			Stopwatch time = new Stopwatch ();
			time.Start ();
			var reads = this.createSequenceProducer (this.Filename);
			var res = SNPCaller.CallSNPs (reads );
			time.Stop ();
			if (this.Verbose) {
				Output.WriteLine (OutputLevel.Verbose);
				Output.WriteLine (OutputLevel.Verbose, "Call Pile-up SNPs: {0}", time.Elapsed);
			}
			return res;

		}

		/// <summary>
		/// Writes the contigs to the file.
		/// </summary>
		/// <param name="assembly">IDeNovoAssembly parameter is the result of running De Novo Assembly on a set of two or more sequences. </param>
		protected void writeContigs (IDeNovoAssembly assembly)
		{
			if (assembly.AssembledSequences.Count == 0) {
				Output.WriteLine (OutputLevel.Results, "No sequences assembled.");
				return;
			}
			ensureContigNames (assembly.AssembledSequences);
            if (!string.IsNullOrEmpty (this.OutputFile)) {
				using (FastAFormatter formatter = new FastAFormatter (this.OutputFile)) {
					formatter.AutoFlush = true;
					foreach (ISequence seq in assembly.AssembledSequences) {
						formatter.Write (seq);
					}
				}
				Output.WriteLine (OutputLevel.Information, "Wrote {0} sequences to {1}", assembly.AssembledSequences.Count, this.OutputFile);
			} else {
				Output.WriteLine (OutputLevel.Information, "Assembled Sequence Results: {0} sequences", assembly.AssembledSequences.Count);
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
			string filename = Path.GetFileNameWithoutExtension (this.OutputFile);
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
		private IEnumerable<CompactSAMSequence> createSequenceProducer (string fileName, DepthOfCoverageGraphMaker coveragePlotter = null)
		{
            if (MakeDepthOfCoveragePlot && !Helper.IsBAM(fileName))
            {
                MakeDepthOfCoveragePlot = false;
                Output.WriteLine(OutputLevel.Error, "Warning: No coverage plots can be made without an input BAM File");
            }

			var parser = new BAMSequenceParser (Filename);
			IEnumerable<CompactSAMSequence> sequences;
			if (ChromosomeName != string.Empty)
            {
				parser.ChromosomeToGet = ChromosomeName;
				sequences = parser.Parse();
            }            
            else
            {
                sequences = parser.Parse();
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
