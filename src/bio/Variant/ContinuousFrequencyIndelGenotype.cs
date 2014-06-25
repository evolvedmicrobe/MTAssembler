using System;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Bio.Variant
{
	public class ContinuousFrequencyIndelGenotype : ContinuousFrequencyGenotype
	{
		public ContinuousFrequencyIndelGenotype ()
		{

		}

		#region implemented abstract members of ContinuousFrequencyGenotype

		public override List<string> GetGenotypesPresent ()
		{
			throw new NotImplementedException ();
		}

		public override List<double> GetGenotypeFrequencies ()
		{
			throw new NotImplementedException ();
		}
		public override string GetMostFrequentGenotype(){
			throw new NotImplementedException ();
		}

		public override double GetHighestFrequency (){
			throw new NotImplementedException ();
		}

		#endregion
	}
}

