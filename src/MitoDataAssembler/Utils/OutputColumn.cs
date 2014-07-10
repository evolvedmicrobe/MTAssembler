using System;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly;
using Bio.Algorithms.Assembly.Padena;
using Bio.Algorithms.Assembly.Graph;
using System.Diagnostics;

namespace MitoDataAssembler
{
	public sealed class OutputColumn<T>
	{
		public readonly string Name;
		public readonly Func<T, object> outFunc;
		public OutputColumn(string name, Func<T, object> outputFunction)
		{
			this.Name = name;
			this.outFunc = outputFunction;
		}

		public static string SafeGet(Func<T, object> func, T GC)
		{
			try
			{
				object o = func(GC).ToString();
				if (o is double)
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


	}
}

