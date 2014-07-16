using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace HaploGrepSharp
{
    public static class CONSTANTS
    {
        /// <summary>
        /// Return true if it passes the filter
        /// </summary>
        /// <param name="poly"></param>
        /// <returns></returns>
        public delegate bool PolymorphismFilterDelegate(Polymorphism poly);
        /// <summary>
        /// A delegate to filter all the polymorphisms
        /// </summary>
        /// <param name="polys"></param>
        /// <returns></returns>
        public delegate IEnumerable<Polymorphism> FilterPolymorphismDelegate(IEnumerable<Polymorphism> polys);

		private static string addFileNameToEndOfDirectory(string fname)
		{
			string rFile = System.Reflection.Assembly.GetExecutingAssembly().Location;
			char splitChar = '\\';
			if (Bio.CrossPlatform.Environment.RunningInMono) {
				splitChar = '/';
			}
			var tempSplit=rFile.Split(splitChar);
			rFile=String.Join("/",tempSplit.Take(tempSplit.Length-1))+splitChar.ToString()+fname;
			return rFile;

		}

		//public static string TREE_XML_FILE { get { return addFileNameToEndOfDirectory ("phylotree15.xml"); } }
		//public static string WEIGHT_FILE { get { return addFileNameToEndOfDirectory ("fluctRates15.txt"); } }
        
		/// <summary>
        /// Applies the current best practices filters to a list of polymorphisms, to 
        /// avoid erroneous calls.  Currently removes indels and bad positions
        /// </summary>
        /// <param name="polys"></param>
        /// <returns></returns>
        internal static IEnumerable<Polymorphism> CommonPolymorphismFilter(IEnumerable<Polymorphism> polys)
        {
            return polys.Where(x => !EXCLUDED_POSITIONS.Contains(x.position) && MutationAssigner.MutationIsBasePair(x.mutation));
        }

        /// <summary>
        /// There are several positions that appear to be problematic for assigning haplotypes 
        /// (See the one-note file and the excel sheet in the source directory for more information.)
        /// I am excluding these at present, and they are attached here.
        /// </summary>
        internal static HashSet<int> EXCLUDED_POSITIONS = new HashSet<int>() {
            302,
            308,
            309,
            310,
            314,
            315,
            455,
            573,
            960,
            965,
            1719,
            2232,
            5899,
            8278,
            515,
            516,
            517,
            518,
            519,
            520,
            521,
            522,
            523,
            524,
            525,
            526,
            8270,
            8271,
            8272,
            8273,
            8274,
            8275,
            8276,
            8277,
            8278,
            8279,
            8280,
            8281,
            8282,
            8283,
            8284,
            8285,
            8286,
            8287,
            15944,
            15945,
            16182,
            16183,
            16193,
            16519};
    }
}
