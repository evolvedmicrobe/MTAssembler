using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using HaploGrepSharp.NewSearchMethods;

namespace HaploGrepSharp.Annotation
{
    public class HaploTypeMatchAnnotation : VariantAnnotation
    {
        HaplotypeComparison bestMatch;
        Dictionary<int, Polymorphism> validPolys;
        string headerLine;
        public HaploTypeMatchAnnotation(HaplotypeComparison bestMatch, Dictionary<int,Polymorphism> validPolys)
        {
            this.bestMatch = bestMatch;
            this.validPolys=validPolys;
            headerLine = String.Join(FIELD_DELIM, new[] { "POSITION", "REF", "SAMPLE", "IN_HAPLOTYPE" });
        }
        public override string GetAnnotation(Polymorphism poly)
        {
            int pos = poly.position;
            string reference, alt, inBestHaplotype;
            alt = "ERROR";
            //TODO: Move this logic in to polymorphism class
            if (poly.mutation == Mutations.INS)
            {
                reference = "-";
            }
            else
            {
                reference = Polymorphism.rCRS[poly.position - 1].ToString();
                alt = poly.InsertedPolys;
            }
            if (poly.mutation == Mutations.DEL)
            {
                alt = "-";
            }
            else if (MutationAssigner.MutationIsBasePair(poly.mutation))
            {
                alt = MutationAssigner.getBase(poly.mutation);
            }
            if (!validPolys.ContainsKey(poly.position))
            {
                inBestHaplotype = "MutationExcludedFromSearch";
            }
            else
            {
                inBestHaplotype = bestMatch.MatchingPolymorphisms.Contains(poly).ToString();
            }
            return String.Join(FIELD_DELIM, new[] { pos.ToString(), reference, alt, inBestHaplotype });
        }
        public override string GetHeaderLine()
        {
            return headerLine;
        }
    }
}
