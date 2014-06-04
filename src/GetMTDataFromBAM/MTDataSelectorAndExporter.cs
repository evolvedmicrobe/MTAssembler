using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio;
using Bio.IO;
using Bio.IO.FastA;
using System.IO;
using Bio.IO.BAM;

namespace GetMTDataFromBAM
{
    public class MTDataSelectorAndExporter
    {
        public bool Verbose;
        public bool CRSAlignedOnly;
        public string OutputFile;
        public string Filename;
        private double pfractionToOutput = .1;
        public string FractionToOutput
        {
            set
            {
                pfractionToOutput = Convert.ToDouble(value);
            }
        }
        public const string BAM_FILE_SUFFIX=".bam";
        public const string DEFAULT_EXPORT_SUFFIX = "_MT.fa";
        public MTDataSelectorAndExporter()
        {
        }
        public void OutputMTReads()
        {
            if(String.IsNullOrEmpty(Filename))
            {
                throw new ArgumentNullException("No input file specified");
            }
            if (!Filename.EndsWith(BAM_FILE_SUFFIX))
            {
                throw new ArgumentNullException("Input file must be a .BAM file");
            }
            if (string.IsNullOrEmpty(OutputFile))
            {
                OutputFile = Filename.Remove(Filename.Length - BAM_FILE_SUFFIX.Length) + DEFAULT_EXPORT_SUFFIX;
            }
            IEnumerable<ISequence> mtReads;
            if(CRSAlignedOnly)
                mtReads= MitoDataGrabber.OutputMitoReadsFromBamFileAlignedToCRSOnly(Filename,pfractionToOutput);
            else
                mtReads=MitoDataGrabber.OutputMitoReadsFromBamFile(Filename);
            
            FastAFormatter fao = new FastAFormatter(OutputFile);
            long count = 0;
            foreach (var seq in mtReads)
            {
                count++;
                fao.Write(seq);
            }
            fao.Close();
            FileInfo FI = new FileInfo(OutputFile);
            
            Console.WriteLine("Wrote "+ count.ToString() + " reads to output file.");
            Console.WriteLine("Of Size: "+GetMTDataFromBAM.Program.FormatMemorySize(FI.Length));
        }
    }
}
