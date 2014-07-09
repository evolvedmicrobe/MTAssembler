using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MitoDataAssembler.IndelCaller
{
    public class IndelVariant
    {
        public int ReferencePosition;
        public List<string> AltAlleles;
        public List<IndelPathCollection.IndelData> OriginalPaths;


    }
}
