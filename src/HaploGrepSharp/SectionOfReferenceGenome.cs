using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio;

namespace HaploGrepSharp
{
    public class SectionOfReferenceGenome
    {
        public bool WrapsAround
        {
            get { return End < Start; }
        }
        public int Start, End;
        public Sequence Seq;

    }
}
