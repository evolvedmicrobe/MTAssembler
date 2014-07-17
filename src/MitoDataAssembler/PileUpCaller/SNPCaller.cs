using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;
using Bio.Variant;
using Bio.Util;
using HaploGrepSharp;
using HaploGrepSharp.NewSearchMethods;

namespace MitoDataAssembler
{
    public class SNPCaller
    {
		/// <summary>
		/// Call SNPs from the sorted list of sequences using the pile-up method.
		/// </summary>
		/// <returns>The SN ps.</returns>
		/// <param name="sequences">Sequences.</param>
		public static SNPCallerReport CallSNPs(IEnumerable<CompactSAMSequence> sequences)
		{
			// Get a pile up and convert it to genotypes
			var pileups = PileUpProducer.CreatePileupFromReads(sequences);
			var genotypes = ContinuousGenotypeCaller.CallContinuousGenotypes (pileups).ToList();

            // Filter down to a usable set
            var usable = genotypes.Where(x => x.ResultType == GenotypeCallResult.GenotypeCalled && x.OriginalPosition.HasValue).ToList();
            if (usable.Count == 0)
            {
				return new SNPCallerReport(AlgorithmResult.Failed);
            }

            // Get median coverage at sites
            var data_counts =usable.Select(x => x.TotalObservedBases).ToList();
            data_counts.Sort();
            var median = data_counts[data_counts.Count / 2];

            //now create a cut-off for required coverage as the square root of the median.
            var cut_off = Math.Sqrt(median);            

			//Get a list of genotypes, and if a simple SNP, make a polymorphism if it doesn't match
			//the reference
			var genotypedPositions = new HashSet<int> ();
			List<Polymorphism> polys = new List<Polymorphism> ();
			foreach (var geno in usable) {
				if (geno.TotalObservedBases >= cut_off) {
					genotypedPositions.Add (geno.OriginalPosition.Value);
					var org_bp = ReferenceGenome.GetReferenceBaseAt_rCRSPosition (geno.OriginalPosition.Value);
					var cur_bp = geno.GetMostFrequentGenotype ();
					if (org_bp != cur_bp[0]) {
						var poly = new Polymorphism (geno.OriginalPosition.Value, MutationAssigner.getBase (cur_bp));
						polys.Add (poly);
					}
				}
			}
			//Now assign haplotype
			HaplotypeSearcher hts = new HaplotypeSearcher ();
			PolymorphismFilter pf = new PolymorphismFilter (p => p.IsSNP && genotypedPositions.Contains (p.Position));
			var simpSample = new SimpleSample ("Pileup", polys, pf);
			var hap_report = hts.GetHaplotypeReport (simpSample);
			return new SNPCallerReport (genotypes, hap_report);
		}
    }
}
