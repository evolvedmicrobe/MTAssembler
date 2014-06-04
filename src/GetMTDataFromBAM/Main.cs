using System;
using SamUtil;
using Bio;
using Bio.IO;
using Bio.IO.BAM;
using System.IO;
using Bio.IO.FastA;
using System.Linq;
using System.Threading.Tasks;
using Bio.Util.ArgumentParser;
using System.Diagnostics;


namespace GetMTDataFromBAM
{
	class Program
	{
        private const double KB = 1024;
        private const double MB = KB * KB;
        private const double GB = MB * KB;
        static void ReportCount(string name, int Count)
        {
            Console.WriteLine(name+" sequences found:\t"+Count.ToString());
        }
		//const string HOME = @"/home/unix/ndelaney/MTData/";
        const string HOME = @"/broad/moothalab/sandbox/ndelaney/";
		const string Data = @"/humgen/1kg/DCC/ftp/data/";
        static void Main(string[] args)
        {
            try
            {
                if ((args == null) || (args.Length < 1))
                {
                    Console.WriteLine("You must provide input arguments");
                }
                else
                {
                    CommandLineArguments clp = new CommandLineArguments();
                    AddParameters(clp);
                    MTDataSelectorAndExporter exporter = new MTDataSelectorAndExporter();
                    clp.Parse(args, exporter);
                    exporter.OutputMTReads();
                }
            }
            catch (Exception ex)
            {
                CatchInnerException(ex);
            }

            Process p = Process.GetCurrentProcess();
            Console.WriteLine("Peak Working Set: ", FormatMemorySize(p.PeakWorkingSet64));
            Console.WriteLine("Total Processor Time: " , p.TotalProcessorTime);
            Debug.Print("Peak Virtual Memory Size 64: ", FormatMemorySize(p.PeakVirtualMemorySize64));
            Debug.Print("Peak Paged Memory Size 64: ", FormatMemorySize(p.PeakPagedMemorySize64));
            
            
           

        }

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
                Console.WriteLine("Error: " + ex.Message);
            }
            else
            {
                CatchInnerException(ex.InnerException);
            }
        }
        private static void AddParameters(CommandLineArguments parser)
        {
            // Add the parameters to be parsed
            parser.Parameter(ArgumentType.Optional, "OutputFile", ArgumentValueType.String, "o", "Output file");
            parser.Parameter(ArgumentType.Optional, "Verbose", ArgumentValueType.Bool, "v", "Display verbose logging during processing.");
            parser.Parameter(ArgumentType.Optional, "CRSAlignedOnly", ArgumentValueType.Bool, "crs", "BAM File was only aligned to CRS");
            parser.Parameter(ArgumentType.Optional, "FractionToOutput", ArgumentValueType.OptionalString, "frac", "Probability each read is included");
            parser.Parameter(ArgumentType.DefaultArgument, "Filename", ArgumentValueType.String, "", "Input file of reads");
            
        }
        static void OldMain(string[] args)
        {
            StreamWriter SW = new StreamWriter(HOME + "CountsByDate.csv");
            StreamReader SR = new StreamReader("FileLocations.csv");
            string[] lines = SR.ReadToEnd().Split('\n');
            lines = lines.Skip(1).ToArray();
            //  Parallel.ForEach(lines, line =>
            Console.WriteLine("Starting");
            foreach (string line in lines)
            {
                Console.WriteLine(line);

                try
                {
                    string[] split = line.Split(',');
                    string fname = split[0];
                    string patientid = split[1];
                    string date = split[2];
                    var mtReads = MitoDataGrabber.OutputMitoReadsFromBamFile(fname);
                    FastAFormatter fao = new FastAFormatter(HOME + patientid + ".fa");
                    long count = 0;
                    foreach (var seq in mtReads)
                    {
                        count++;
                        fao.Write(seq);
                    }
                    fao.Close();
                    FileInfo FI = new FileInfo(fname);
                    string size = FI.Length.ToString();
                    lock (SW)
                    {
                        SW.WriteLine(String.Join(",", patientid, count.ToString(), size, date));
                        Console.WriteLine(patientid + " has " + count.ToString() + " reads");
                    }
                    if (args.Length > 2)
                        break;
                }
                catch (Exception thrown)
                { Console.WriteLine(thrown.Message); }
            }
            //);
            SW.Close();
        }

    }
			
           
        //string f = d + "NA21144.mapped.ILLUMINA.bwa.GIH.exome.20111114.bam";
				// BAMParser bp = new BAMParser();
				// var map = bp.ParseRange(f, "MT");
				//BAMFormatter bamf = new BAMFormatter();
				//bamf.CreateIndexFile = true;
				//bamf.Format(map, d + "Test.bam");
				//foreach (DirectoryInfo DI in new DirectoryInfo(Data).GetDirectories())
				//{
					//string newDir=HOME+DI.Name+"/";
					//string oldDir=DI.FullName+"/exome_alignment/";
					//Directory.CreateDirectory(newDir);
					//DirectoryInfo DI2 = new DirectoryInfo(oldDir);
					//if(DI2.Exists && !DI.Name.Contains("chrom") && !DI.Name.Contains("SOLID"))
					//{
						//foreach (FileInfo FI in DI2.GetFiles())
					//{
						//if (FI.Name.EndsWith(".bam"))
						//{
}
