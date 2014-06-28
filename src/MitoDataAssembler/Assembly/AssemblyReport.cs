using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly.Padena;
using System.Reflection;

namespace MitoDataAssembler
{
	public sealed class AssemblyReport : AlgorithmReport
    {
        public static FieldInfo[] outputValues;
        
		static AssemblyReport()
        {
            outputValues = typeof(MitoPaintedAssembler).GetFields().Where(prop => prop.IsDefined(typeof(OutputAttribute), true)).ToArray();
	     }
        public static string CreateHeaderLine()
        {
            return String.Join(",", outputValues.Select(x => x.Name).ToArray());
        }
		/// <summary>
		/// Creates an an empty assembly report if the method was never attempted.
		/// </summary>
		public AssemblyReport() : base(AlgorithmResult.NotAttempted)
		{
			HeaderLineForCSV = CreateHeaderLine ();
			this.DataLineForCSV = String.Join (",", Enumerable.Repeat ("NA", outputValues.Length));
		}

		public AssemblyReport(MitoPaintedAssembler toReportOn, AlgorithmResult result = AlgorithmResult.Success) : base(result)
        {
			HeaderLineForCSV = CreateHeaderLine ();

			this.DataLineForCSV =  String.Join(",", outputValues.Select(x => x.GetValue(toReportOn).ToString()).ToArray());
        }       

    } 

}
