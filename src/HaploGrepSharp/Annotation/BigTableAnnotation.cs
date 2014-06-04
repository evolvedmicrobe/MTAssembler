using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace HaploGrepSharp.Annotation
{
    /// <summary>
    /// This is an annotation file made by (I believe) Steve that for each position lists the 
    /// possible substition, and it's effect, etc.
    /// </summary>
    public class BigTableAnnotation : VariantAnnotation
    {
        string header;
        const int headerOutputColsSkip=4;
        Dictionary<string, BigTableData> positionToData;
        string nullAnnotation;
        /// <summary>
        /// Data from one data line
        /// </summary>
        public class BigTableData
        {
            /// <summary>
            /// Non coding regions have this "snp" annotation
            /// </summary>
            public const string ALL_SITE_SAME="*";
            public int Position;
            public string RefBP;
            public string AltBP;
            public string Data;
            public BigTableData(string line)
            {
                var sp = line.Split('\t');
                Position = Convert.ToInt32(sp[1]);
                var pos2=Convert.ToInt32(sp[2]);
                if(Position!=pos2)
                {
                    throw new HaploGrepException("Annotation made a bad assumption about the file");
                }
                RefBP=sp[3];
                AltBP=sp[4];
                Data=String.Join(FIELD_DELIM,sp.Skip(headerOutputColsSkip).ToArray());                
            }
        }

        protected string makeQueryString(int position,string refBP,string altBP)
        {
            return position.ToString()+refBP+altBP;
        }

        public BigTableAnnotation()
        {
            string data = Properties.Resources.crs_rn_bigtable;
            var lines = data.Split('\n');
            header=String.Join(FIELD_DELIM,lines[0].Split('\t').Skip(headerOutputColsSkip).ToArray());
            var temp=lines.Skip(1).Where(y=>!String.IsNullOrEmpty(y)).Select(x=>new BigTableData(x)).ToList();
            //positionToData=temp.ToDictionary(x=>makeQueryString(x.Position,x.RefBP,x.AltBP));
            positionToData = new Dictionary<string, BigTableData>();
            foreach (var v in temp)
            {
                var key=makeQueryString(v.Position,v.RefBP,v.AltBP);
                //TODO: Merge overlapping annotations
                if(!positionToData.ContainsKey(key)) {
                    positionToData[key] = v;
                }
            }
            nullAnnotation=createNAstring(header.Split(FIELD_DELIM[0]).Length);
        }
        public override string GetHeaderLine()
        {
            return header;
        }
        public override string GetAnnotation(Polymorphism p)
        {
            if (MutationAssigner.MutationIsComplex(p.Mutation))
            {
                return nullAnnotation;
            }
            else
            {
                var site = p.Position;
                var refBP=Polymorphism.getReferenceBaseSingle(site);
                var alt=MutationAssigner.getBase(p.Mutation);
                var key = makeQueryString(site,refBP,alt);
                if (!positionToData.ContainsKey(key))
                {
                    key=makeQueryString(site,BigTableData.ALL_SITE_SAME,BigTableData.ALL_SITE_SAME);
                    if (positionToData.ContainsKey(key))
                    {
                        return positionToData[key].Data;
                    }
                    return nullAnnotation;
                }                
                else
                {
                    return positionToData[key].Data;
                }
            }
        }
    }
}
