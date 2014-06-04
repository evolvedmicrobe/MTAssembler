using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace HaploGrepSharp.Annotation
{
    /// <summary>
    /// Annotates a variant position based on the frequency file Danny compiled from
    /// http://www.mtdb.igp.uu.se/
    /// </summary>
    public class FrequencyFileAnnotation  : VariantAnnotation
    {
        string fname="mtDB.sorted.crs.rn.txt";
        Dictionary<int, double> positionToFrequency=new Dictionary<int,double>();
        public FrequencyFileAnnotation()
        {
            string data = Properties.Resources.mtDB_sorted_crs_rn;
            var lines = data.Split('\n');
            foreach(var l in lines.Skip(1))
            {
                if (String.IsNullOrEmpty(l)) continue;
                var sp = l.Split('\t');
                var pos = Convert.ToInt32(sp[0].Replace("CRS:", ""));
                var freq = Convert.ToDouble(sp[1]);
                positionToFrequency[pos] = freq;
            }
            

        }
        public override string GetHeaderLine()
        {
            return "MTDB_Mut_Freq";
        }
        public override string GetAnnotation(Polymorphism p)
        {
            if (MutationAssigner.MutationIsComplex(p.Mutation))
            {
                return "NA";
            }
            else
            {
                var site = p.Position;
                if (!positionToFrequency.ContainsKey(site))
                {
                    return "NA";
                }
                else
                {
                    return positionToFrequency[site].ToString();
                }
            }
        }
    }
}
