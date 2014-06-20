using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using Bio.IO.BAM;
using Bio.IO.SAM;
using Bio;
using System.Diagnostics;

namespace MitoDataAssembler
{
    /// <summary>
    /// A class designed to create a depth of coverage plot based on a BAM file for the MT DNA.
    /// NOT THREAD SAFE.
    /// </summary>
    public class DepthOfCoverageGraphMaker
    {
        double[] depthCounts = new double[StaticResources.CRS_LENGTH];
        byte gap = DnaAlphabet.Instance.Gap;    
        public DepthOfCoverageGraphMaker()
        {

        }       
        /// <summary>
        /// Unpacks the sequence so that it is aligned to the reference at the given start but ignoring insertions.
        /// Useful for Depth of Coverage;
        /// </summary>
        public void ProcessCountCoverageFromSequence(CompactSAMSequence orgSeq)
        {
			if (orgSeq != null && orgSeq.RName != StaticResources.MT_CHROMOSOME_NAME)
            { return;}

            string CIGAR = orgSeq.CIGAR;
            if (string.IsNullOrWhiteSpace(CIGAR) || CIGAR=="*")
            {
                return;
            }
            int curRef = orgSeq.Pos-1; 
            List<KeyValuePair<char, int>> charsAndPositions = new List<KeyValuePair<char, int>>();
            for (int i = 0; i < CIGAR.Length; i++)
            {
                char ch = CIGAR[i];
                if (Char.IsDigit(ch))
                {
                    continue;
                }
                charsAndPositions.Add(new KeyValuePair<char, int>(ch, i));
            }
            for (int i = 0; i < charsAndPositions.Count; i++)
            {
                int start = 0;
                int end = 0;
                int len = 0;
                char ch = charsAndPositions[i].Key;
                if (i == 0)
                {
                    start = 0;
                }
                else
                {
                    start = charsAndPositions[i - 1].Value + 1;
                }
                end = charsAndPositions[i].Value - start;
                len = int.Parse(CIGAR.Substring(start, end));
                switch (ch)
                {
                    case 'P': //padding (Silent deltions from padded reference)
                    case 'N': //skipped region from reference
                        throw new Exception("Not built to handle clipping yet");                        
                    case 'M': //match or mismatch
                    case '=': //match
                    case 'X': //mismatch
                        for (int k = 0; k < len; k++)
                        {
                            if (curRef >= StaticResources.CRS_LENGTH)
                            {
								Debug.WriteLine ("Seq: " + orgSeq.ID + " is aligned past the MT DNA reference genome");
                                break;
                            }
							depthCounts [curRef] = depthCounts [curRef] + 1.0;; 
                            curRef++;                            
                        }
                        break;
                    case 'I'://insertion to the reference
						break;
                    case 'D'://Deletion from the reference
						curRef += len;
						break;
                    case 'S': //soft clipped
                    case 'H'://had clipped
                        break;
                    default:
                        throw new FormatException("Unexpected SAM Cigar element found " + ch.ToString());
                }
            }
        }
        /// <summary>
        /// Make a PDF depth of coverage plot for the given file and MT DNA name.
        /// Can be called after the reads are all processed.
        /// </summary>
        /// <param name="BAMFileName"></param>
        /// <param name="outputFileName"></param>
        /// <returns></returns>
        public void OutputCoverageGraphAndCSV(string outputFilePrefix)
        {
            //Now to make the plot
            var r = new RInterface();
            List<RInterface.OptionalArgument> args = new List<RInterface.OptionalArgument>() { new RInterface.OptionalArgument("type", "l") };
			r.PlotPDF(Enumerable.Range(0, StaticResources.CRS_LENGTH).Select(x => (double)(x+1)), depthCounts, outputFilePrefix+"_Coverage.pdf", "Depth of Coverage", "mtDNA Position", "Depth",null,args);
            StreamWriter sw = new StreamWriter(outputFilePrefix + "_Coverage.csv");
            sw.WriteLine("Pos,Depth");
            for (int i = 0; i < depthCounts.Length; i++)
            {
                sw.WriteLine((i + 1).ToString() + "," + depthCounts[i].ToString());
            }
            sw.Close();
        }
    }
}
