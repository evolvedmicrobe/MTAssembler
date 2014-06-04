using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HaploGrepSharp.NewSearchMethods
{    public class SimpleSample
    {
        //TODO: Make this static?
        public readonly PolymorphismFilter Filter;
        public string Id;
        public List<Polymorphism> Polymorphisms;
        public double TotalSampleWeight { get; private set; }
        public SimpleSample(string id, List<Polymorphism> polys, PolymorphismFilter filter)
        {
            this.Filter=filter;
            this.Id=id;
            //TODO: See if the filter was already applied?
            this.Polymorphisms=filter.FilterPolys(polys).ToList();
            //now calculate weight
            this.TotalSampleWeight = Polymorphisms.Sum(x => x.getMutationRate());
        }
    }
}
