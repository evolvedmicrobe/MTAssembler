using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio.IO;
using Bio;
using System.Collections;
using System.Collections.Generic;
using Bio.IO.SAM;

namespace MitoDataAssembler
{
    public static class ReadFilter
    {
		static bool FilterDuplicates = true;//FilterDuplicates
        static bool TrimQuality = true;
        static private int _trimEndQuality=22;
        static private int _meanRequiredQuality=22;        
		public static IEnumerable<CompactSAMSequence> FilterReads(IEnumerable<CompactSAMSequence> preFiltered, DepthOfCoverageGraphMaker coverageCounter = null)
        {
            foreach (var toFilter in preFiltered)
            {
				if ( FilterDuplicates && ((toFilter.SAMFlags & SAMFlags.Duplicate) == SAMFlags.Duplicate))
                { continue; }   

				//Process coverage before trimming, as otherwise the CIGARs are off...
				if (coverageCounter != null) {
					coverageCounter.ProcessCountCoverageFromSequence (toFilter);
				}
                
                int[] vals = toFilter.GetQualityScores();
				int lastAcceptableBase = (int)toFilter.Count - 1;
                while (lastAcceptableBase >= 0)
                {
                    if (vals[lastAcceptableBase] >= _trimEndQuality)
                    {
                        break;
                    }
                    lastAcceptableBase--;
                }
                if (lastAcceptableBase > 0)
                {
                    //check mean
                    double mean = vals.Take(lastAcceptableBase + 1).Average();
                    if (mean > _meanRequiredQuality)
                    {
						//only trim if necessary.
						if (lastAcceptableBase < (toFilter.Count-1)) {
							yield return toFilter.GetSubSequence (0, lastAcceptableBase + 1) as CompactSAMSequence;
						} else {
							yield return toFilter;
						}
					 }
                }
                
            }
        }

    }
}
