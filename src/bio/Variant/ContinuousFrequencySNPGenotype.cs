using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Bio.Variant
{
	public class ContinuousFrequencySNPGenotype : ContinuousFrequencyGenotype
	{


		public BasePairFrequencies Frequencies;

		public ContinuousFrequencySNPGenotype(BasePairFrequencies frequencies) {
			ResultType = GenotypeCallResult.GenotypeCalled;
			Frequencies = Frequencies;
		}

		public ContinuousFrequencySNPGenotype(GenotypeCallResult res) {
			ResultType = res;
		}


		#region implemented abstract members of ContinuousFrequencyGenotype
		static List<string> basesInOrder = Enumerable.Range(0,4).
				Select(x => BaseAndQuality.Get_DNA_MappingForIndex(x).
				ToString()).
				ToList();
		public override List<string> GetGenotypesPresent ()
		{
			return basesInOrder;
		}

		public override List<double> GetGenotypeFrequencies ()
		{
			return Frequencies.Frequencies.ToList ();
		}

		#endregion
	}
}

