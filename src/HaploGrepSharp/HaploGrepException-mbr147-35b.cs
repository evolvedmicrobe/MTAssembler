using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace HaploGrep
{
	/// <summary>
	/// The exception for haplogrep classes
	/// </summary>
    public class HaploGrepException : Exception
    {
        public HaploGrepException(string message) : base(message) { }
        public HaploGrepException(string message, Exception innerException) : base(message, innerException) { }
    }
}
