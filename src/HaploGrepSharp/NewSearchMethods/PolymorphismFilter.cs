using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HaploGrepSharp.NewSearchMethods
{
    public class PolymorphismFilter
    {
        CONSTANTS.PolymorphismFilterDelegate additionalFilter;
        public PolymorphismFilter() {
            additionalFilter = null;
        }
        public PolymorphismFilter(CONSTANTS.PolymorphismFilterDelegate additionalFilter)
        {
            this.additionalFilter = additionalFilter;          
        }
        public IEnumerable<Polymorphism> FilterPolys(IEnumerable<Polymorphism> polys)
        {
            if (additionalFilter == null){
                return CONSTANTS.CommonPolymorphismFilter(polys);
            }
            else
                return CONSTANTS.CommonPolymorphismFilter(polys.Where(z => additionalFilter(z)));
        }
    }
}
