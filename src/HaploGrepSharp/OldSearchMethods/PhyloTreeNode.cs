using System.Collections.Generic;

namespace HaploGrepSharp
{

    public class PhyloTreeNode
    {
        internal Haplogroup haplogroup = null;
        internal List<Polymorphism> expectedPolys = new List<Polymorphism>();
        internal List<Polymorphism> foundPolys = new List<Polymorphism>();
        internal List<Polymorphism> notInRangePolys = new List<Polymorphism>();
        internal List<Polymorphism> correctedBackmutation = new List<Polymorphism>();
        public PhyloTreeNode(Haplogroup haplogroup)
        {
            this.haplogroup = haplogroup;

        }
        public virtual void addExpectedPoly(Polymorphism currentPoly)
        {
            this.expectedPolys.Add(currentPoly);

        }
        public virtual void addCorrectPoly(Polymorphism currentPoly)
        {
            this.foundPolys.Add(currentPoly);

        }
        public virtual void addNotInRangePoly(Polymorphism currentPoly)
        {
            this.notInRangePolys.Add(currentPoly);

        }
        public virtual void removeExpectedPoly(Polymorphism currentPoly)
        {
            this.expectedPolys.Remove(currentPoly);

        }
        public virtual void removeCorrectPoly(Polymorphism currentPoly)
        {
            this.foundPolys.Remove(currentPoly);

        }

        public virtual Haplogroup Haplogroup
        {
            get
            {
                return this.haplogroup;

            }
        }

        public virtual List<Polymorphism> ExpectedPolys
        {
            get
            {
                return this.expectedPolys;

            }
        }

        public virtual List<Polymorphism> FoundPolys
        {
            get
            {
                return this.foundPolys;

            }
        }

        public virtual List<Polymorphism> NotInRangePolys
        {
            get
            {
                return this.notInRangePolys;

            }
        }

        public virtual List<Polymorphism> CorrectedBackmutation
        {
            get
            {
                return this.correctedBackmutation;

            }
        }
        public virtual void addCorrectedBackmutation(Polymorphism poly)
        {
            this.correctedBackmutation.Add(poly);

        }
    }
}