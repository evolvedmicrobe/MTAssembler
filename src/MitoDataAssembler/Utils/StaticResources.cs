using System;

namespace MitoDataAssembler
{
	public static class StaticResources
	{
        /// <summary>
        /// Out of date, and not currently used.
        /// </summary>
		public const string ASSEMBLE_HELP = @"Usage:	 MTAssembler.exe assemble [options] -ref:<ReferenceGenome(fasta)> <Read File> 		
Description: Assemble reads into contigs (No matepair information required.) 

Typically Used Parameters:
-k:<int>           Length of k-mer.
-o:<string>        Output file
-chr:<string>	   Typicallly 'MT' or the name of chromosome to select
				   from.
-v                 Display verbose logging during processing.
-p:<string>        Output a collection of diagnostic information
                   files with the given prefix.


Infrequently Used Parameters:
-h                 Help
-q                 Display minimal output during processing.
-d:<int>           Threshold for removing 
                   dangling ends in graph.
-r:<int>    	   Length Threshold for 
                   removing redundant paths in graph.
-e:<int>           Threshold for eroding low coverage ends.
-i                 Bool to do erosion or not. 
-a                 Whether to estimate kmer length.
-c:<int>           Threshold used for removing 
                   low-coverage contigs.
";

		
		public const string ContinuePrompt = "Are you sure you want to continue (Y/N)?";
		
		public const string PeakPagedMemorySize64 = "Peak memory in the virtual memory paging file used: {0}";
			
		public const string PeakVirtualMemorySize64 = "Peak virtual memory used: {0}";
			
		public const string PeakWorkingSet64 = "Peak physical memory used: {0}";
			
		public const string ReferenceFile = "Reference File not in correct format";

		public const string TotalProcessorTime = "Total CPU time taken: {0}";
			
		public const string UnknownCommand = "Unknown command: {0}";
        
        /// <summary>
        /// The expected name of the mtDNA chromosome in the BAM File.
        /// </summary>
        public const string MT_CHROMOSOME_NAME = "MT";

        /// <summary>
        /// The length of the CRS chromosome.
        /// </summary>
        public const int CRS_LENGTH = 16568;

        /// <summary>
        /// Deletions less than this value are considered indels, above this value are large deletions.
        /// </summary>
        public const int SIZE_DIF_BETWEEN_LARGE_AND_SMALL_DELETION = 150;
	}
}

