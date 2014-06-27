using System;

namespace MitoDataAssembler
{
	/// <summary>
	/// A general class that returns results from different stages of the algorithm, to
	/// be merged later on.
	/// </summary>
	public class AlgorithmReport
	{
        public AlgorithmResult Result { get; protected set; }

		/// <summary>
		/// To simplify results, one string is returned from each BAM processed that can be added to a CSV.
		/// This string is the concatenation of each analysis subset, and this field gives the header lines for it.
		/// </summary>
		public string HeaderLineForCSV { get; protected set;}

		/// <summary>
		/// The data line for the CSV report string.
		/// </summary>
		public string DataLineForCSV{ get; protected set;}
		public AlgorithmReport (AlgorithmResult result)
		{
            Result = result;
		}
	}

    public enum AlgorithmResult
    {
        Failed, Success
    }
}

