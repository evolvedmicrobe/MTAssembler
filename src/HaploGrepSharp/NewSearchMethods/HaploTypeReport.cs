using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HaploGrepSharp.NewSearchMethods
{
	public class ReportColumn
	{
		public delegate string GetReportValue (HaploTypeReport rpt);

		public static string ReportDelimiter = ",";
		public readonly string ColumnName;
		GetReportValue valueGetter;

		public ReportColumn (string name, GetReportValue func)
		{
			this.ColumnName = name;
			this.valueGetter = func;
		}

		public string GetValue (HaploTypeReport report)
		{
			try {
				return valueGetter (report);
			} catch {
				return "NA";
			}
		}
	}

	public class HaploTypeReport
	{
		public HaplotypeComparison BestHit;
		public int NumberOfEquallyGoodBestHits;
		public HaplotypeComparison SecondBestHit;
		public SimpleSample Sample;

		public HaploTypeReport (List<HaplotypeComparison> Comparisons, SimpleSample sample)
		{
			Sample = sample;
			var temp = Comparisons.Where (x => x.node.haplogroup.id == "U5b").First ();

			Comparisons.Sort ((x, y) => -x.Rank.CompareTo (y.Rank));
			BestHit = Comparisons [0];
			NumberOfEquallyGoodBestHits = Comparisons.TakeWhile (x => x.Rank == BestHit.Rank).Count ();
			SecondBestHit = Comparisons [1];
		}

		public IEnumerable<string> GetRowReportLines ()
		{
			yield return "Best Haplotype = " + BestHit.node.haplogroup.id;
			yield return "Best Score = " + BestHit.Rank.ToString ();
			yield return "Total Haplotypes with Best Score = " + NumberOfEquallyGoodBestHits.ToString ();
			yield return "Second Best Haplotype = " + SecondBestHit.node.haplogroup.id;
			yield return "Second Best Score = " + SecondBestHit.Rank.ToString ();
			yield return "Matched Polymorphisms = " + BestHit.NumberOfMatchingPolymorphisms.ToString ();
			yield return "Polymorphisms missing from Haplotype = " + BestHit.NumberOfPolymorphismsMissingFromHaplotype.ToString ();
			yield return "Polymorphisms missing from Genotype = " + BestHit.NumberOfPolymorphismsMissingFromGenotype.ToString ();
		}

		public string GetColumnReportLine ()
		{
			return String.Join (ReportColumn.ReportDelimiter, ReportValues.Select (x => x.GetValue (this)).ToArray ());
		}

		public static string GetColumnReportHeaderLine (string prefix = "")
		{
			return String.Join (ReportColumn.ReportDelimiter, ReportValues.Select (x => prefix+x.ColumnName).ToArray ());
		}

		static List<ReportColumn> ReportValues = new List<ReportColumn> () { 
			new ReportColumn ("SampleName", x => x.Sample.Id),
			new ReportColumn ("BestHaplotype", x => x.BestHit.node.haplogroup.id),
			new ReportColumn ("BestScore", x => x.BestHit.Rank.ToString ()),
			new ReportColumn ("EquallyBestHits", x => x.NumberOfEquallyGoodBestHits.ToString ()),
			new ReportColumn ("SecondBestHaplotype", x => x.SecondBestHit.node.haplogroup.id),
			new ReportColumn ("SecondBestScore", x => x.SecondBestHit.Rank.ToString ()),
			new ReportColumn ("TotalSampleDifferencesFromRCRS", x => x.Sample.Polymorphisms.Count.ToString ()),
			new ReportColumn ("TotalSampleWeight", x => x.Sample.TotalSampleWeight.ToString ()),
			new ReportColumn ("PolymorphismsInBestHaplotype", x => x.BestHit.NumberOfPolymorhpismsInHaplotype.ToString ()),
			new ReportColumn ("PolymorphismsMissingFromBestHaplotype", x => x.BestHit.NumberOfPolymorphismsMissingFromHaplotype.ToString ()),
			new ReportColumn ("PolymorphismsMissingFromGenotype", x => x.BestHit.NumberOfPolymorphismsMissingFromGenotype.ToString ())
		};
	}
}
