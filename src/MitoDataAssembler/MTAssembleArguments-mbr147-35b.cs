﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using Bio;
using Bio.Algorithms.Assembly;
using Bio.Algorithms.Assembly.Padena;
using Bio.IO;
using Bio.IO.FastA;
using Bio.Util;

namespace MitoDataAssembler
{
    /// <summary>
    /// Class for Assemble options.
    /// </summary>
    public class MTAssembleArguments
    {
        #region Public Fields
        /// <summary>
        /// Length of k-mer.
        /// </summary>
        public int KmerLength = 19;

		public string ReferenceGenome;

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

        /// <summary>
        /// Output file.
        /// </summary>
        public string OutputFile = string.Empty;

        /// <summary>
        /// If a subset of the BAM is used this can be specified here.
        /// </summary>
        public string ChromosomeName = string.Empty;

        /// <summary>
        /// Prefix of the report output prefix file
        /// </summary>
        public string ReportOutputPrefix = "AssemblyReport_";

        /// <summary>
        /// Display verbose logging during processing.
        /// </summary>
        public bool Verbose
        {
            get { return _pVerbose; }
            set
            {
                if (value)
                    Output.TraceLevel = OutputLevel.Information | OutputLevel.Verbose;
                else
                    Output.TraceLevel = OutputLevel.Information;
                _pVerbose = value;
                //TODO: Make this an option again
                _pVerbose = true;

            }
        }
        private bool _pVerbose=true;
      
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
        /// It assembles the sequences.
        /// </summary>
        public virtual void AssembleSequences()
        {
            AssembleSequencesReturningString();
            
        }
        /// <summary>
        /// Assembles the sequences and returns the string that can be placed in a CSV output report.
        /// </summary>
        /// <returns></returns>
        public string AssembleSequencesReturningString()
        {
            string toReturn = "No String Set";
            TimeSpan algorithmSpan = new TimeSpan();
            Stopwatch runAlgorithm = new Stopwatch();
            FileInfo refFileinfo = new FileInfo(this.Filename);
            long refFileLength = refFileinfo.Length;

            runAlgorithm.Restart();
            IEnumerable<ISequence> reads = this.ParseFile();
            runAlgorithm.Stop();
            algorithmSpan = algorithmSpan.Add(runAlgorithm.Elapsed);

            Output.WriteLine(OutputLevel.Information, StaticResources.AssemblyStarting);

            if (this.Verbose)
            {
                Output.WriteLine(OutputLevel.Verbose);
                Output.WriteLine(OutputLevel.Verbose, "Processed read file: {0}", Path.GetFullPath(this.Filename));
                Output.WriteLine(OutputLevel.Verbose, "   Read/Processing time: {0}", runAlgorithm.Elapsed);
                Output.WriteLine(OutputLevel.Verbose, "   File Size           : {0}", refFileLength);
                Output.WriteLine(OutputLevel.Verbose, "   k-mer Length        : {0}", this.KmerLength);
            }

            using (MitoPaintedAssembler assembler = new MitoPaintedAssembler())
            {
                assembler.AllowErosion = true;
                //assembler.ReferenceGenomeFile = ReferenceGenome;
                assembler.DiagnosticFileOutputPrefix = DiagnosticFilePrefix;
				Console.WriteLine("Prefix is: "+assembler.DiagnosticFileOutputPrefix);
				Console.WriteLine ("Diagnostic Information On: "+assembler.OutputDiagnosticInformation.ToString ());

                assembler.StatusChanged += this.AssemblerStatusChanged;
                assembler.AllowErosion = this.AllowErosion;
                assembler.AllowKmerLengthEstimation = this.AllowKmerLengthEstimation;

                if (ContigCoverageThreshold != -1)
                {
                    assembler.AllowLowCoverageContigRemoval = true;
                    assembler.ContigCoverageThreshold = ContigCoverageThreshold;
                }
                assembler.DanglingLinksThreshold = this.DangleThreshold;
                assembler.ErosionThreshold = this.ErosionThreshold;
                if (!this.AllowKmerLengthEstimation)
                {
                    assembler.KmerLength = this.KmerLength;
                }

                assembler.RedundantPathLengthThreshold = this.RedundantPathLengthThreshold;
                runAlgorithm.Restart();
                IDeNovoAssembly assembly = assembler.Assemble(reads);
                 runAlgorithm.Stop();
                algorithmSpan = algorithmSpan.Add(runAlgorithm.Elapsed);
                if (this.Verbose)
                {
                    Output.WriteLine(OutputLevel.Verbose);
                    Output.WriteLine(OutputLevel.Verbose, "Compute time: {0}", runAlgorithm.Elapsed);
                }

                runAlgorithm.Restart();
                this.WriteContigs(assembly);
                runAlgorithm.Stop();
                toReturn = assembler.GetReportLine();
                if (assembler.OutputDiagnosticInformation)
                {
                    var outFile = new StreamWriter(ReportOutputPrefix + DiagnosticFilePrefix + ".csv");
                    outFile.WriteLine(MitoDataAssembler.AssemblyOutput.CreateHeaderLine());
                    outFile.WriteLine(toReturn);
                    outFile.Close();

 
                }
                algorithmSpan = algorithmSpan.Add(runAlgorithm.Elapsed);
                if (this.Verbose)
                {
                    Output.WriteLine(OutputLevel.Verbose);
                    Output.WriteLine(OutputLevel.Verbose, "Write contigs time: {0}", runAlgorithm.Elapsed);
                    Output.WriteLine(OutputLevel.Verbose, "Total runtime: {0}", algorithmSpan);
                }
            }
            return toReturn;
        }
        #endregion

        #region Protected Members

        /// <summary>
        /// It Writes the contigs to the file.
        /// </summary>
        /// <param name="assembly">IDeNovoAssembly parameter is the result of running De Novo Assembly on a set of two or more sequences. </param>
        protected void WriteContigs(IDeNovoAssembly assembly)
        {
            if (assembly.AssembledSequences.Count == 0)
            {
                Output.WriteLine(OutputLevel.Results, "No sequences assembled.");
                return;
            }

            EnsureContigNames(assembly.AssembledSequences);

            if (!string.IsNullOrEmpty(this.OutputFile))
            {
                using (FastAFormatter formatter = new FastAFormatter(this.OutputFile))
                {
                    formatter.AutoFlush = true;

                    foreach (ISequence seq in assembly.AssembledSequences)
                    {
                        formatter.Write(seq);
                    }
                }
                Output.WriteLine(OutputLevel.Information, "Wrote {0} sequences to {1}", assembly.AssembledSequences.Count, this.OutputFile);
            }
            else
            {
                Output.WriteLine(OutputLevel.Information, "Assembled Sequence Results: {0} sequences", assembly.AssembledSequences.Count);
                using (FastAFormatter formatter = new FastAFormatter())
                {
                    formatter.Open(new StreamWriter(Console.OpenStandardOutput()));
                    formatter.MaxSymbolsAllowedPerLine = decideOutputWidth() ;
                    formatter.AutoFlush = true;
                    foreach (ISequence seq in assembly.AssembledSequences)
                    {
                        formatter.Write(seq);
                    }
                }
            }
        }
		private int decideOutputWidth()
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
        private void EnsureContigNames(IList<ISequence> sequences)
        {
            for (int index = 0; index < sequences.Count; index++)
            {
                ISequence inputSequence = sequences[index];
				if (string.IsNullOrEmpty (inputSequence.ID))
					inputSequence.ID = GenerateSequenceId (index + 1);
				else
					inputSequence.ID = GenerateSequenceId (index + 1) + inputSequence.ID;
            }
        }

        /// <summary>
        /// Generates a sequence Id using the output filename, or first input file.
        /// </summary>
        /// <param name="counter">Sequence counter</param>
        /// <returns>Auto-generated sequence id</returns>
        private string GenerateSequenceId(int counter)
        {
            string filename = Path.GetFileNameWithoutExtension(this.OutputFile);
            if (string.IsNullOrEmpty(filename))
                filename = Path.GetFileNameWithoutExtension(this.Filename);
            filename = filename.Replace(" ", "");
            Debug.Assert(!string.IsNullOrEmpty(filename));
            return string.Format(CultureInfo.InvariantCulture, "{0}_{1}", filename, counter);
        }

        /// <summary>
        /// Parses the File.
        /// </summary>
        protected IEnumerable<ISequence> ParseFile()
        {
            return ParseFile(this.Filename);
        }

        /// <summary>
        /// Helper method to parse the given filename into a set
        /// of ISequence elements. 
        /// </summary>
        /// <param name="fileName">Filename to load data from</param>
        /// <returns>Enumerable set of ISequence elements</returns>
        internal IEnumerable<ISequence> ParseFile(string fileName)
        {
            ISequenceParser parser = SequenceParsers.FindParserByFileName(fileName);
            if (parser == null)
                parser = new FastAParser(fileName);
            //else if (parser is Bio.IO.FastQ.FastQZippedParser)
            //{
            //    Console.WriteLine("Filtering FastQ File");
            //    var filter = new Bio.Filters.FastQZippedQualityTrimmedParser(fileName);
            //    return  ReadFilter.FilterReads(filter.Parse());
            //}
            if(parser is Bio.IO.BAM.BAMSequenceParser &&  ChromosomeName!=string.Empty)
            {
                var b = parser as Bio.IO.BAM.BAMSequenceParser;
                b.ChromosomeToGet = ChromosomeName;
                
            }
            
            return ReadFilter.FilterReads(parser.Parse());
        }
        /// <summary>
        /// Method to handle status changed event.
        /// </summary>
        protected void AssemblerStatusChanged(object sender, StatusChangedEventArgs statusEventArgs)
        {
            if (Verbose)
                Output.WriteLine(OutputLevel.Verbose,  statusEventArgs.StatusMessage);
            else if (!Quiet)
            {
                if (statusEventArgs.StatusMessage.StartsWith("Step", StringComparison.OrdinalIgnoreCase)
                    && statusEventArgs.StatusMessage.Contains("Start"))
                {
                    int pos = statusEventArgs.StatusMessage.IndexOf(" - ", StringComparison.OrdinalIgnoreCase);
                    Output.WriteLine(OutputLevel.Information, statusEventArgs.StatusMessage.Substring(0,pos));
                }
            }
        }     
        #endregion
    }
}