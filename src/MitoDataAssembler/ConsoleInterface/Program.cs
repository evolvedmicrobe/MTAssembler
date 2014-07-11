using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using Bio;
using Bio.IO;
using Bio.Util.ArgumentParser;
using Bio.Util;

namespace MitoDataAssembler
{
    public class Program
    {
        private const double KB = 1024;
        private const double MB = KB * KB;
        private const double GB = MB * KB;


    /// <summary>
    /// Main function of program class.
    /// </summary>
    /// <param name="args">Arguments to the main function.</param>
    public static void Main(string[] args)
    {
        args=args.Select(x=>x.Replace("%SKYDRIVE%",Environment.GetEnvironmentVariable("SKYDRIVE"))).ToArray();
        //Takes two arguments, the first is the reference genome, the second is the reads file.				
        Output.WriteLine(OutputLevel.Required, "\nThe Mitochondrial Assembler ");
		if (Bio.CrossPlatform.Environment.RunningInMono) {
				Console.WriteLine ("Running in the Mono runtime");
		}
            DateTime compTime = RetrieveLinkerTimestamp();
			Console.WriteLine("Compiled on: " +compTime.ToLongDateString()+" "+compTime.ToLongTimeString());
#if !DEBUG
            try
            {
#endif
				ProcessData(args);
#if !DEBUG
            }

            catch (Exception ex)
            {
                CatchInnerException(ex);
            }
#endif
                Process p = Process.GetCurrentProcess();
                Console.WriteLine("\nResource Report:");
		        Console.WriteLine(StaticResources.PeakWorkingSet64, FormatMemorySize(p.PeakWorkingSet64));
                Console.WriteLine(StaticResources.TotalProcessorTime, p.TotalProcessorTime);
                Console.WriteLine(StaticResources.PeakVirtualMemorySize64, FormatMemorySize(p.PeakVirtualMemorySize64));
                Console.WriteLine(StaticResources.PeakPagedMemorySize64, FormatMemorySize(p.PeakPagedMemorySize64));
			}

    #region Private Methods
		/// <summary>
		/// Formats the specified memory in bytes to appropriate string.
		/// for example, 
		///  if the value is less than one KB then it returns a string representing memory in bytes.
		///  if the value is less than one MB then it returns a string representing memory in KB.
		///  if the value is less than one GB then it returns a string representing memory in MB.
		///  else it returns memory in GB.
		/// </summary>
		/// <param name="value">value in bytes</param>
		public static string FormatMemorySize(long value)
		{
			string result = null;
			if (value > GB)
			{
				result = (Math.Round(value / GB, 2)).ToString() + " GB";
			}
			else if (value > MB)
			{
				result = (Math.Round(value / MB, 2)).ToString() + " MB";
			}
			else if (value > KB)
			{
				result = (Math.Round(value / KB, 2)).ToString() + " KB";
			}
			else
			{
				result = value.ToString() + " Bytes";
			}

			return result;
		}

		/// <summary>
		/// Catches Inner Exception Messages.
		/// </summary>
		/// <param name="ex">Exception</param>
		private static void CatchInnerException(Exception ex)
		{
			if (ex.InnerException == null || string.IsNullOrEmpty(ex.InnerException.Message))
			{
				Output.WriteLine(OutputLevel.Error, "Error: " + ex.Message);
                Output.WriteLine(OutputLevel.Error, ex.StackTrace);
		        throw ex;
			}
			else
			{
				CatchInnerException(ex.InnerException);
			}
		}
	
		/// <summary>
		/// Main function to do all the processing.
		/// </summary>
		/// <param name="args">Arguments to Assemble.</param>
		private static void ProcessData(string[] args)
		{
			MTAssembleArguments options = new MTAssembleArguments();
			CommandLineArguments parser = new CommandLineArguments();
			AddAssembleParameters(parser);
			if (args.Length < 1 
			|| args[0].Equals("Help", StringComparison.InvariantCultureIgnoreCase)
			|| args[0].Equals("/h", StringComparison.CurrentCultureIgnoreCase)
			|| args[0].Equals("/help", StringComparison.CurrentCultureIgnoreCase)
			|| args[0].Equals("-h", StringComparison.CurrentCultureIgnoreCase)
			|| args[0].Equals("-help", StringComparison.CurrentCultureIgnoreCase))
			{
				//Output.WriteLine(OutputLevel.Required, StaticResources.ASSEMBLE_HELP);
				PrintHelp (parser);
			}
			else
			{
				try
				{
					parser.Parse(args, options);
				}
				catch (ArgumentParserException ex)
				{
					Output.WriteLine(OutputLevel.Error, ex.Message);
					//Output.WriteLine(OutputLevel.Required, StaticResources.ASSEMBLE_HELP);
					Environment.Exit(-1);
				}
				if (options.Help)
				{
					PrintHelp (parser);
				}
				else
				{
					if (options.Verbose)
						Output.TraceLevel = OutputLevel.Information | OutputLevel.Verbose;
					else if (!options.Quiet)
						Output.TraceLevel = OutputLevel.Information;
					options.ProcessMTDNA();
				}
			}				
		}
        private static void AddAssembleParameters(CommandLineArguments parser)
		{
			// Add the parameters to be parsed
			parser.Parameter (ArgumentType.DefaultArgument, "Filename", ArgumentValueType.String, "", "Input BAM file with reads");
			//parser.Parameter (ArgumentType.Optional, "OutputFile", ArgumentValueType.String, "o", "Output file");
			parser.Parameter (ArgumentType.Optional, "Help", ArgumentValueType.Bool, "h", "Display help");        

			parser.Parameter (ArgumentType.Optional, "Quiet", ArgumentValueType.Bool, "q", "Display minimal output during processing.");
			parser.Parameter (ArgumentType.Optional, "KmerLength", ArgumentValueType.OptionalInt, "k", "Length of k-mer");
			parser.Parameter (ArgumentType.Optional, "DangleThreshold", ArgumentValueType.OptionalInt, "d", "Threshold for removing dangling ends in graph");
			parser.Parameter (ArgumentType.Optional, "RedundantPathLengthThreshold", ArgumentValueType.OptionalInt, "r", "Length Threshold for removing redundant paths in graph");
			parser.Parameter (ArgumentType.Optional, "ErosionThreshold", ArgumentValueType.OptionalInt, "e", "Threshold for eroding low coverage ends");
			parser.Parameter (ArgumentType.Optional, "AllowErosion", ArgumentValueType.Bool, "i", "Bool to do erosion or not.");
			parser.Parameter (ArgumentType.Optional, "AllowKmerLengthEstimation", ArgumentValueType.Bool, "a", "Whether to estimate kmer length.");
			parser.Parameter (ArgumentType.Optional, "ContigCoverageThreshold", ArgumentValueType.Int, "c", "Threshold used for removing low-coverage contigs.");
			parser.Parameter (ArgumentType.Optional, "Verbose", ArgumentValueType.Bool, "v", "Display verbose logging during processing.");
			parser.Parameter (ArgumentType.Optional, "Skip_DepthOfCoveragePlot", ArgumentValueType.Bool, "skip_dp", "Skip depth of coverage plot");
			//parser.Parameter (ArgumentType.Required, "ReferenceGenome", ArgumentValueType.String, "ref","Reference Genome File (Fasta");
			parser.Parameter (ArgumentType.Optional, "DiagnosticFilePrefix", ArgumentValueType.String, "o", "Prefix to append to all diagnostic files.");
			parser.Parameter (ArgumentType.Optional, "ChromosomeName", ArgumentValueType.String, "chr", "Name of mtDNA Chromsome (typicallly rCRS or MT).");
			parser.Parameter (ArgumentType.Optional, "Skip_Pileup_Calling", ArgumentValueType.Bool, "skip_pileup", "Skip calling SNPs and haplotypes using a columnwise pile-up in addition to the de novo assembly");
			parser.Parameter (ArgumentType.Optional, "Skip_EM_Frequency_Estimates", ArgumentValueType.Bool, "skip_em", "Estimate SNP frequencies by simple read counts.");
			parser.Parameter (ArgumentType.Optional, "Skip_Assembly_Step", ArgumentValueType.Bool, "skip_asm", "Skip the de-novo assembly and indel finding step");
			parser.Parameter (ArgumentType.Optional, "Skip_Peak_Finder", ArgumentValueType.Bool, "skip_pf", "Skip the peak finding step");
            parser.Parameter(ArgumentType.Optional, "ForceSqrtThreshold", ArgumentValueType.Bool, "force_sqrt", "Always uses Sqrt(Median Coverage) as cut-off");
            parser.Parameter(ArgumentType.Optional, "MinimumNodeCount", ArgumentValueType.Int, "min_node_count", "Erode graph nodes with coverage <= this value.");
			parser.Parameter(ArgumentType.Optional, "OutputIntermediateGraphSteps", ArgumentValueType.Bool, "show_graphs", "Output graphs and coverage histograms at all assembly steps.");
            parser.Parameter(ArgumentType.Optional, "NoContigOutput", ArgumentValueType.Bool, "no_contigs", "Skip output of assembly contigs.");

            
        }

        private static DateTime RetrieveLinkerTimestamp()
            {
                string filePath = System.Reflection.Assembly.GetCallingAssembly().Location;
                const int c_PeHeaderOffset = 60;
                const int c_LinkerTimestampOffset = 8;
                byte[] b = new byte[2048];
                System.IO.Stream s = null;

                try
                {
                    s = new System.IO.FileStream(filePath, System.IO.FileMode.Open, System.IO.FileAccess.Read);
                    s.Read(b, 0, 2048);
                }
                finally
                {
                    if (s != null)
                    {
                        s.Close();
                    }
                }

                int i = System.BitConverter.ToInt32(b, c_PeHeaderOffset);
                int secondsSince1970 = System.BitConverter.ToInt32(b, i + c_LinkerTimestampOffset);
                DateTime dt = new DateTime(1970, 1, 1, 0, 0, 0);

                dt = dt.AddSeconds(secondsSince1970);
                dt = dt.AddHours(TimeZone.CurrentTimeZone.GetUtcOffset(dt).Hours);
                return dt;
            }
		private static void PrintHelp(CommandLineArguments parser)
		{
			Output.WriteLine(OutputLevel.Required,@"
MTAssembler.exe assemble [options] <Input BAM File> 
Description: Assemble reads into contigs (No matepair information required.)");
			Output.Write (OutputLevel.Required, parser.PrintArguments ());
		}
			#endregion
	}
}
