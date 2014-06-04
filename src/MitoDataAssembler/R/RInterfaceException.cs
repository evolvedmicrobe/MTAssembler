using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MitoDataAssembler.R
{
    public class RInterfaceException :Exception
    {
        public RInterfaceException(string msg) : base(msg) { }

    }
}
