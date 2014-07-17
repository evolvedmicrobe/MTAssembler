using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using HaploGrepSharp;
using HaploGrepSharp.TreeUtilities;

namespace HaploGrepSharp.NewSearchMethods
{
	/// <summary>
	/// Class represents a comparison between a particular haplotype and 
	/// a sample
	/// </summary>
	public class HaplotypeComparison
	{
		/// <summary>
		/// The node in phylotree to compare a sample against
		/// </summary>
		public PhyloTreeNodev2 node;

		/// <summary>
		/// The rank by haplogrep, where higher is better
		/// </summary>
		public double Rank { get; private set; }

		/// <summary>
		/// The weight of all polymorphisms in the node that pass the filter, formerly usedWeightPolys
		/// </summary>
		double haplotypeWeightPoly;
		/// <summary>
		/// The weight of all matching mutations, previously correctWeightPolys
		/// </summary>
		double matchingWeightPoly;
		/// <summary>
		/// All polymorphisms shared between the two.
		/// </summary>
		public List<Polymorphism> MatchingPolymorphisms;

		public int NumberOfMatchingPolymorphisms { get { return MatchingPolymorphisms.Count; } }

		public ushort NumberOfPolymorphismsMissingFromHaplotype;
		public ushort NumberOfPolymorphismsMissingFromGenotype;
		public ushort NumberOfPolymorhpismsInHaplotype;

		/// <summary>
		/// Create a new class that compares a sample to a known haplotype
		/// </summary>
		public HaplotypeComparison (PhyloTreeNodev2 node, SimpleSample sample)
		{
			this.node = node;
			var polysInNode = sample.Filter.FilterPolys (node.Mutations).ToList ();
			NumberOfPolymorhpismsInHaplotype = (ushort)polysInNode.Count;
			haplotypeWeightPoly = polysInNode.Sum (x => x.getMutationRate ());
			MatchingPolymorphisms = polysInNode.Where (z => sample.Polymorphisms.Contains (z)).ToList ();
            NumberOfPolymorphismsMissingFromGenotype = (ushort)(polysInNode.Count - MatchingPolymorphisms.Count); 
            NumberOfPolymorphismsMissingFromHaplotype= (ushort)(sample.Polymorphisms.Count - MatchingPolymorphisms.Count);
			matchingWeightPoly = MatchingPolymorphisms.Sum (x => x.getMutationRate ());
			var CorrectInHaplotypeRatio = haplotypeWeightPoly == 0 ? 1.0 : (matchingWeightPoly / haplotypeWeightPoly);
			var CorrectInSampleRatio = sample.TotalSampleWeight == 0 ? 1.0 : (matchingWeightPoly / sample.TotalSampleWeight);
			Rank = .5 * CorrectInHaplotypeRatio + .5 * CorrectInSampleRatio;
		}
	}
}
