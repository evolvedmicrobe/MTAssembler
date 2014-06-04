using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using HaploGrepSharp.NewSearchMethods;
using HaploGrepSharp.TreeUtilities;
using Bio;
using HaploGrepSharp;
using HaploGrepSharp.Annotation;

namespace HaploGrepSharp.NewSearchMethods
{
	/// <summary>
	/// The basic class for finding haplotypes given a complete mtDNA sequence
	/// can be extended later to account for SNP Chips
	/// </summary>
	public class HaplotypeSearcher
	{
		List<PolymorphismFilter> filtersApplied = new List<PolymorphismFilter> ();
		PhyloTreeNodev2 treeRoot;

		public HaplotypeSearcher ()
		{
			treeRoot = TreeLoader.LoadTree (CONSTANTS.TREE_XML_FILE);
			var skipped = HaploSearchManager.SetPhylogeneticWeights (CONSTANTS.WEIGHT_FILE);

		}

		public HaploTypeReport GetHaplotypeReport (Sequence toBuildReportFor, List<string> dataLines, string id = "Sample")
		{
			var delts = ReferenceGenome.GetDeltaAlignments (toBuildReportFor).SelectMany (x => x).ToList ();
			if (delts.Count != 1) {
                return null;
				throw new Exception ("Final assembly had no or multiple delta alignments with the reference, whereas only 1 is expected");
			}
			var delt = delts [0];
			var aln = delt.ConvertDeltaToSequences ();
			dataLines.Add ("REF_OFFSET =" + delt.FirstSequenceStart.ToString ());
			dataLines.Add ("ASSEMBLY_OFFSET =" + delt.SecondSequenceStart.ToString ());
			//dataLines.Add("ALN_SCORE=" + aln.Score.ToString());
			var refseq = aln.FirstSequence as Sequence;
			dataLines.Add ("REF_SEQUENCE  = " + refseq.ConvertToString ());
			var qseq = aln.SecondSequence as Sequence;
			dataLines.Add ("QUERY_SEQUENCE= " + qseq.ConvertToString ());            
            
			//now get all polymorphisms, and sort haplotypes
			//TODO: Don't call nucmer twice
			var AllPolys = GenomeToHaploGrepConverter.FindPolymorphisms (toBuildReportFor);
			PolymorphismFilter pf = new PolymorphismFilter ();
			var sample = new SimpleSample (id, AllPolys, pf);
			var comparisons = treeRoot.GetAllChildren ().Select (x => new HaplotypeComparison (x, sample)).ToList ();
			var report = new HaploTypeReport (comparisons, sample);
			dataLines.AddRange (report.GetRowReportLines ());           
			var passedFilters = sample.Polymorphisms.ToDictionary (x => x.position);

			//noq get list of differences
            dataLines.Add("#Data report below, (Note in rCRS positions, accounting for 'N' at position 3107 and in 1 based index)");
            List<VariantAnnotation> reporters=new List<VariantAnnotation>() {new HaploTypeMatchAnnotation(report.BestHit,passedFilters),
                                                                              new FrequencyFileAnnotation(),new BigTableAnnotation()};
            var headerLine=String.Join(VariantAnnotation.FIELD_DELIM,reporters.Select(x=>x.GetHeaderLine()).ToArray());
			dataLines.Add (headerLine);
            var infoLines=AllPolys.Select(x=> String.Join(VariantAnnotation.FIELD_DELIM,reporters.Select(y=>y.GetAnnotation(x)).ToArray()));
            dataLines.AddRange(infoLines);
			return report;
		}

		/// <summary>
		/// TODO: Should verify that this only happens once, perhaps also that the samples don't contain things that wouldn't match this, 
		/// </summary>
		/// <param name="filterToApply">Filter to apply.</param>
		public void ApplyFilterToAllSNPs (PolymorphismFilter filterToApply)
		{
			foreach (var root in treeRoot.GetAllChildren ()) {
				root.Mutations.FilterSites (filterToApply);
			}
			filtersApplied.Add (filterToApply);
		}

		public HaploTypeReport GetHaplotypeReport (SimpleSample sample)
		{
			//verify filtering is the same
			foreach (var filt in filtersApplied) {
				if (sample.Polymorphisms.Count () != filt.FilterPolys (sample.Polymorphisms).Count ()) {
					throw new HaploGrepException ("It appears a different filter was used when constructing the tree as opposed to " +
					"constructing the samples");
				}
			}
			var comparisons = treeRoot.GetAllChildren ().Select (x => new HaplotypeComparison (x, sample)).ToList ();
			var report = new HaploTypeReport (comparisons, sample);
			return report;

		}
	}
}

