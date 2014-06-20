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
        static bool Filter1024 = true;//FilterDuplicates
        static bool TrimQuality = true;
        static private int _trimEndQuality=22;
        static private int _meanRequiredQuality=22;        
        public static IEnumerable<ISequence> FilterReads(IEnumerable<ISequence> preFiltered, DepthOfCoverageGraphMaker coverageCounter = null)
        {
            foreach (var toFilter in preFiltered)
            {
                var samRead = toFilter as CompactSAMSequence;
				if ( samRead !=null)
                {
                    if ((samRead.SAMFlags & 1024) == 1024)
                    { continue; }   
					if (coverageCounter != null) {
						coverageCounter.ProcessCountCoverageFromSequence (samRead);
					}
                }
                var qs = toFilter as QualitativeSequence;
                if (qs != null)
                {
                    int[] vals = qs.GetQualityScores();
                    int lastAcceptableBase = (int)qs.Count - 1;
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
                            yield return qs.GetSubSequence(0, lastAcceptableBase + 1) as QualitativeSequence;
                        }
                    }
                }
                else
                {
                    yield return toFilter;
                }
            }
        }

    }
}
