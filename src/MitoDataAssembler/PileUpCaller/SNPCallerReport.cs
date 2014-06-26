using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Variant;
using HaploGrepSharp.NewSearchMethods;

namespace MitoDataAssembler
{
    public class SNPCallerReport : AlgorithmReport
    {
		public List<ContinuousFrequencyGenotype> Genotypes { get; private set; }

		public HaploTypeReport HaplotypeInformation { get; private set; }

		public SNPCallerReport(List<ContinuousFrequencyGenotype> genos, HaploTypeReport hap_report)
		{
			this.Genotypes = genos;
			this.HaplotypeInformation = hap_report;

			HeaderLineForCSV = HaploTypeReport.GetColumnReportHeaderLine ("snp_");
			this.DataLineForCSV = hap_report.GetColumnReportLine ();
		}



    }
}
