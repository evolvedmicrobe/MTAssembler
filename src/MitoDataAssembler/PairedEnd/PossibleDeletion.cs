using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MitoDataAssembler.PairedEnd
{
    public class PossibleDeletion
    {
        public int Start, End;
    }
    public class PossibleDeletionCollection : List<PossibleDeletion>
    {
        public PossibleDeletionCollection() : base() { }
        public override string ToString()
        {
            if (this.Count == 0) { return "NONE"; }
            else
            {
                StringBuilder sb=new StringBuilder();
                int i=1;
                foreach (var pd in this)
                {
                    sb.Append(i.ToString() + ":" + pd.Start.ToString() + "-" + pd.End.ToString()+"; ");
                }
                return sb.ToString();
            }
        }

    }
}
