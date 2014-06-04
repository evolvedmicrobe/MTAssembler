using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace HaploGrepSharp.Annotation
{
    /// <summary>
    /// Abstract class for adding annotations to finished assemblies
    /// </summary>
    public abstract class VariantAnnotation
    {
        public static string FIELD_DELIM="\t";
        public abstract string GetHeaderLine();
        public abstract string GetAnnotation(Polymorphism p);
        /// <summary>
        /// Creates a list of NA's to return if the annotation is irrelevant.
        /// </summary>
        /// <param name="numOfNA"></param>
        /// <returns></returns>
        protected string createNAstring(int numOfNA)
        {
            return String.Join(FIELD_DELIM, Enumerable.Range(0, numOfNA).Select(x => "NA").ToArray());
        }
    }
}
