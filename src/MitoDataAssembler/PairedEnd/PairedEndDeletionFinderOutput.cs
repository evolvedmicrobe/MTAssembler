﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MitoDataAssembler.PairedEnd;
using System.Reflection;
using Bio.Algorithms.Assembly.Padena;

namespace MitoDataAssembler.PairedEnd
{
	public static class PairedEndDeletionFinderOutput
	{
		public static PropertyInfo[] outputValues;

		public static int TotalOutputFields {
			get { return outputValues.Length; }
		}

		static PairedEndDeletionFinderOutput ()
		{
			outputValues = typeof(PairedEndPeakFinder).GetProperties ().Where (prop => prop.IsDefined (typeof(OutputAttribute), true)).ToArray ();
		}

		public static string CreateHeaderLine ()
		{
			return String.Join (",", outputValues.Select (x => x.Name).ToArray ());
		}

		public static string GetReportValues (PairedEndPeakFinder toReportOn)
		{
               
			return String.Join (",", outputValues.Select (x => GetValueSafe (toReportOn, x)).ToArray ());
		}

		public static string GetValueSafe (PairedEndPeakFinder toReportOn, PropertyInfo FI)
		{
			try {
				return FI.GetValue (toReportOn).ToString ();
			} catch {
				return "NULL";
			}
		}
	}
}
