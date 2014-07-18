using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using Bio;
using Bio.IO.BAM;
using Bio.Util;

namespace MitoDataAssembler
{
	public class BAMNuclearChromosomeReadGenerator
	{
		
        /// <summary>
        /// Create a sequence enumerator that filters the reads and adds them to the depth of coverage counter
        /// if necessary.
        /// </summary>
        /// <param name="bamFileName">Filename to load data from</param>
        /// <returns>Enumerable set of ISequence elements</returns>
        public static IEnumerable<CompactSAMSequence> GetNuclearAndMitochondrialReads(string bamFileName)
        {

            var parser = new BAMSequenceParser(bamFileName);
            var header = parser.GetFileHeader();
            foreach (var chr in Regions.Split('\n'))
            {
                var chr2 = chr.Trim();
                var split = chr2.Split('\t');
                var locs = split.Skip(1).Select(x =>Convert.ToInt32(x)).ToList();
                var name = split[0];
                if (!header.ReferenceSequences.Any( z => z.Name == name))
                {
                    name = "chr"+name;
                    if(!header.ReferenceSequences.Any( z => z.Name == name))
                    {
                        continue;
                    }
                }

                foreach (var v in parser.ParseRangeAsEnumerableSequences(bamFileName, name, locs[0], locs[1]))
                {
                    yield return v;
                }  
            }
        }
        /// <summary>
        /// List of regions where mtDNA reads align:
        /// 
        /// http://evolvedmicrobe.com/blogs/?p=234
        /// </summary>
        public static string Regions =
@"1	564461	570304
1	77436909	77436955
11	10530238	10530313
11	10530438	10530564
11	10530721	10531312
11	10531456	10531548
11	10531696	10531768
13	110076479	110076572
13	110076643	110076723
17	22020693	22020787
17	22020838	22020965
17	51183348	51183442
17	51183553	51183666
2	88124412	88124531
2	88124643	88124803
2	149639293	149639421
3	96336531	96336658
3	96336751	96336847
5	8621096	8621132
5	79946379	79946454
5	79946664	79947124
5	79947246	79947657
5	134260858	134261004
5	134261122	134261204
5	134262289	134262368
5	134262469	134262625
5	134262718	134262835
5	134262922	134263058
5	134263625	134263700
5	134263911	134263986
MT	1	16569
X	125864879	125864917";
	}
}

