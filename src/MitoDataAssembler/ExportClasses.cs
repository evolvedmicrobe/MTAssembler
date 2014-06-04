using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MitoDataAssembler
{
    public class Export
    {
        public sealed class OutputColumn
        {
            public readonly string Name;
            public Func<MitoPaintedAssembler, object> outFunc;
            public OutputColumn(string Name, Func<MitoPaintedAssembler, object> OutputFunction)
            {
                this.Name = Name;
                this.outFunc = OutputFunction;
            }
        }
        public static string SafeGet(Func<MitoPaintedAssembler, object> func, MitoPaintedAssembler ass)
        {
            try
            {
                object o= func(ass).ToString();
                if(o is double)
                {
                    double n = (double)o;
                    return n.ToString("n5");
                }
                else
                {
                    return o.ToString();
                }
            }
            catch
            {
                return "NA";
            }
        }
        public static List<OutputColumn> OutputColumnCollection = new List<OutputColumn>() {
                //new OutputColumn("Name", (x) => x.DataSetName),
                //new OutputColumn("Doubling Time(Hrs)", (x) => x.GrowthRate.DoublingTime),              
                };
       
    }
   
}
