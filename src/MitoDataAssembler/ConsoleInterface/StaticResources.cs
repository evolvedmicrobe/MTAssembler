using System;

namespace MitoDataAssembler
{
	public static class StaticResources
	{
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

		

		public const string AssemblyStarting = "Beginning assembly processing.";
				
				public const string BadKmerLength=@"WARNING: kmerLength is set to an even value ({0}) which can lead to infinite loops when building the graphs.
			It is recommended that you use odd numbers to ensure palindrome sequences do not lead to infinite recursion.";
				
				public const string ContinuePrompt
="Are you sure you want to continue (Y/N)?";
				

				public const string PeakPagedMemorySize64
="Peak memory in the virtual memory paging file used: {0}";
			
				public const string PeakVirtualMemorySize64
="Peak virtual memory used: {0}";
			
				public const string PeakWorkingSet64
="Peak physical memory used: {0}";
			
				public const string ReferenceFile
="Reference File not in correct format";
				

				public const string TotalProcessorTime
="Total CPU time taken: {0}";
			
				public const string UnknownCommand="Unknown command: {0}";

	}
}

