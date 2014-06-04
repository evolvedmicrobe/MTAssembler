using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace HaploGrepSharp.Annotation
{
    /// <summary>
    /// Based on a file Sarah put together of the same name as the resource
    /// </summary>
    public class ConfirmedDiseaseAssociation : VariantAnnotation
    {
        Dictionary<int, string> positionToFrequency;
        string noAnnotationResponse;
        string[] headerCols = new[] { "Nuc_change", "AA_change", "Gene", "DiseaseType" };
        public ConfirmedDiseaseAssociation()
        {
            noAnnotationResponse=String.Join(FIELD_DELIM,headerCols.Select(x=>"NA").ToArray());
            string data = Properties.Resources.mitomap_dz_mutations_110126_filtered_110916_confirmed_120126;
            var lines = data.Split('\n');
            foreach(var l in lines.Skip(1))
            {
                var sp = l.Split('\t');
                var pos = Convert.ToInt32(sp[0].Replace("CRS:", ""));
                var annotation=String.Join(FIELD_DELIM,sp.Skip(1).ToArray());
                positionToFrequency[pos] = annotation;
            }
            

        }
        public override string GetHeaderLine()
        {
            return String.Join(FIELD_DELIM,headerCols);
        }
        public override string GetAnnotation(Polymorphism p)
        {
                var site = p.Position;
                if (!positionToFrequency.ContainsKey(site))
                {
                    return noAnnotationResponse;
                }
                else
                {
                    return positionToFrequency[site];
                }
        }
    }
}
