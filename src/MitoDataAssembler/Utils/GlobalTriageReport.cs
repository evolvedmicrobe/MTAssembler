using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MitoDataAssembler.Utils
{
    /// <summary>
    /// This class is for the universal "bubbling up" of messages when errors are thrown in some algorithmic step, but we do 
    /// not want to bring down the entire program because of it.  Usually the erroneous data will be ignored and the algorithm will
    /// continue but it can be reported here to prevent problems.
    /// </summary>
    public static class GlobalTriageReport
    {
        public class ErrorReport
        {
            public string ErrorLocation;
            public string ErrorDescription;
            public ErrorReport(string step, string description) { 
                
            }
        }
        static List<ErrorReport> Errors = new List<ErrorReport>();
        /// <summary>
        /// Add an internal error to the global reportion mechanism, 
        /// so that it is output at the end.
        /// </summary>
        /// <param name="step"></param>
        /// <param name="description"></param>
        public static void ReportError(string step, string description)
        {
            lock (Errors)
            {
                Errors.Add(new ErrorReport(step, description));
            }
        }
    }
}
