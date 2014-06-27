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

		//Actually is one plus this number.
		static private int _minReadLength = 18;
		static private int _trimEndQuality = 22;
		static private int _meanRequiredQuality = 22; 

		public static IEnumerable<CompactSAMSequence> FilterReads(IEnumerable<CompactSAMSequence> preFiltered, DepthOfCoverageGraphMaker coverageCounter = null)
        {
            foreach (var toFilter in preFiltered)
            {
				if ( FilterDuplicates && ((toFilter.SAMFlags & SAMFlags.Duplicate) == SAMFlags.Duplicate))
                { continue; }   

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
				if (lastAcceptableBase > _minReadLength)
                {
                    //check mean
                    double mean = vals.Take(lastAcceptableBase + 1).Average();
                    if (mean > _meanRequiredQuality)
                    {
						//only trim if necessary.
						CompactSAMSequence toReturn;
						if (lastAcceptableBase < (toFilter.Count-1)) {
							toFilter.TrimSequence (lastAcceptableBase + 1);
						}
						//Process coverage before returning the read
						if (coverageCounter != null) {
							coverageCounter.ProcessCountCoverageFromSequence (toFilter);
						}
						yield return toFilter;
					 }
                }

                
            }
        }

    }
}
