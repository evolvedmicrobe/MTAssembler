using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Variant;
using HaploGrepSharp.NewSearchMethods;
using System.IO;

namespace MitoDataAssembler
{
    public class SNPCallerReport : AlgorithmReport
    {
		public List<ContinuousFrequencySNPGenotype> Genotypes { get; private set; }

		public HaploTypeReport HaplotypeInformation { get; private set; }

		public SNPCallerReport(List<ContinuousFrequencySNPGenotype> genos, HaploTypeReport hap_report) :base(AlgorithmResult.Success)
		{
			this.Genotypes = genos;
			this.HaplotypeInformation = hap_report;

			HeaderLineForCSV = HaploTypeReport.GetColumnReportHeaderLine ("snp_");
			DataLineForCSV = hap_report.GetColumnReportLine ();
		}

		public SNPCallerReport(AlgorithmResult result) : base(result) {
			HeaderLineForCSV = HaploTypeReport.GetColumnReportHeaderLine ("snp_");
			DataLineForCSV = String.Join(",",Enumerable.Repeat("NA",HeaderLineForCSV.Count(x=>x==',')-1));            
        }

		public void OutputHaploReport(string prefix) {
			var sw = new StreamWriter(prefix +"SNP_Haplotypes.csv");
			sw.WriteLine (HeaderLineForCSV);
			sw.WriteLine (DataLineForCSV);
			sw.Close ();
		}



    }
}
