using System;
using System.Collections.Generic;
using System.Globalization;

namespace Bio.IO.SAM
{
	public static class CigarUtils
	{
		public static List<CigarElement> GetCigarElements(string cigar)
		{
			var elements = new List<CigarElement> (7);
			int start = 0;
			for (int i = 1; i < cigar.Length; i++) {
				char ch = cigar [i];
				if (!Char.IsDigit (ch)) {
					var str_len = i - start;
					var length = int.Parse (cigar.Substring (start, str_len), CultureInfo.InvariantCulture);
					elements.Add (new CigarElement (ch, length));
					start = i++;//increment i as the next element will be a number
				}
			}
		}
		/// <summary>
		/// These are the CIGAR elements "MDNX=" can change the alignment length on the reference
		/// this check if it is one of these elements.
		/// </summary>
		/// <returns><c>true</c>, if elementis_ MDN x_ equal was cigared, <c>false</c> otherwise.</returns>
		/// <param name="element">Element.</param>
		public static bool CigarElementis_MDNX_Equal(char element)
		{
			//"MDNX=";
			return element == 'M' || element == 'D' || element == 'N' ||
			element == 'X' || element == '=';
		}
	}
}

