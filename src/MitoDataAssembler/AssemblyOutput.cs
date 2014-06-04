using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly.Padena;
using System.Reflection;

namespace MitoDataAssembler
{
    public sealed class AssemblyOutput
    {
        public static FieldInfo[] outputValues;
        static AssemblyOutput()
        {
            outputValues = typeof(MitoPaintedAssembler).GetFields().Where(prop => prop.IsDefined(typeof(OutputAttribute), true)).ToArray();

        }
        public static string CreateHeaderLine()
        {
            return String.Join(",", outputValues.Select(x => x.Name).ToArray());
        }
        public static string GetReportValues(MitoPaintedAssembler toReportOn)
        {
            return String.Join(",", outputValues.Select(x => x.GetValue(toReportOn).ToString()).ToArray());
        }       

    } 

}
