using System;
using System.Collections;
using System.Collections.Generic;
namespace HaploGrepSharp
{
    /// <summary>
    /// A class that contains the results of checking for a polymorphism along one branch of the tree.
    /// This seems to be made largely by the HaploSearchClass
    /// </summary>
    public class SearchResult : IComparable<SearchResult>
    {
        private Haplogroup haplogroup;
        private List<Polymorphism> allCheckedPolys = new List<Polymorphism>();
        private List<Polymorphism> correctPolys = new List<Polymorphism>();
        private List<Polymorphism> unusedPolys = new List<Polymorphism>();
        private List<Polymorphism> unusedPolysNotInRange = new List<Polymorphism>();
        private List<Polymorphism> correctedBackmutations = new List<Polymorphism>();
        private List<Polymorphism> missingPolysOutOfRange = new List<Polymorphism>();
        private PhyloTreePath usedPath = new PhyloTreePath();
        private Sample usedPolysInSample = null;
        private double usedWeightPolys = 0.0D;
        private double correctWeightPolys = 0.0D;
        private double expectedWeightPolys = 0.0D;
        private double missingWeightPolys = 0.0D;

        public SearchResult(string haplogroup,  Sample polysInTestSample)
        {
            this.haplogroup = new Haplogroup(haplogroup);
            this.usedPolysInSample = polysInTestSample;
            this.unusedPolys.AddRange(polysInTestSample.Polymorphismn);
            foreach (Polymorphism currentPoly in polysInTestSample.Polymorphismn) {
                this.usedWeightPolys += currentPoly.getMutationRate();
            }
        }
        public SearchResult(string newHaplogroup, SearchResult resultToCopy) {
            this.haplogroup = new Haplogroup(newHaplogroup);
            this.usedPolysInSample = resultToCopy.usedPolysInSample;
            this.allCheckedPolys.AddRange(resultToCopy.allCheckedPolys);
            this.correctPolys.AddRange(resultToCopy.correctPolys);
            this.unusedPolys.AddRange(resultToCopy.unusedPolys);
            this.correctedBackmutations.AddRange(resultToCopy.correctedBackmutations);
            this.unusedPolysNotInRange.AddRange(resultToCopy.unusedPolysNotInRange);
            this.missingPolysOutOfRange.AddRange(resultToCopy.missingPolysOutOfRange);
            this.usedPath = new PhyloTreePath(resultToCopy.usedPath);
            this.usedWeightPolys = resultToCopy.usedWeightPolys;
            this.correctWeightPolys = resultToCopy.correctWeightPolys;
            this.expectedWeightPolys = resultToCopy.expectedWeightPolys;
            this.missingWeightPolys = resultToCopy.missingWeightPolys;
        }
        public Haplogroup Haplogroup   {
            get  {
                return this.haplogroup;
            }
        }
        public List<Polymorphism> CheckedPolys {
            get {
                return this.allCheckedPolys;
            }
        }
        public double Rank {
            get {
                return CorrectPolyInTestSampleRatio * 0.5D + CorrectPolyInHaplogroupRatio * 0.5D;
            }
        }
        public double CorrectPolyInTestSampleRatio
        {
            get
            {
                return this.correctWeightPolys / this.usedWeightPolys;
            }
        }
        public double CorrectPolyInHaplogroupRatio
        {
            get {
                if (this.expectedWeightPolys != 0.0D) {
                    return this.correctWeightPolys / this.expectedWeightPolys;
                } return 1.0D;
            }
        }
        public  List<Polymorphism> CorrectPolys
        {
            get
            {
                return this.correctPolys;
            }
        }
        public  List<Polymorphism> MissingPolysOutOfRange
        {
            get
            {
                return this.missingPolysOutOfRange;
            }
        }
        public  List<Polymorphism> UnusedPolys
        {
            get
            {
                return this.unusedPolys;
            }
        }
        public  List<Polymorphism> UnusedPolysNotInRange
        {
            get
            {
                return this.unusedPolysNotInRange;
            }
        }
        
        public  Sample Sample
        {
            get
            {
                return this.usedPolysInSample;
            }
        }
        public  double UsedWeightPolys
        {
            get
            {
                return this.usedWeightPolys;
            }
        }
        public  double CorrectWeightPolys
        {
            get
            {
                return this.correctWeightPolys;
            }
        }
        public  double ExpectedWeightPolys
        {
            get
            {
                return this.expectedWeightPolys;
            }
        }
        public  void addCorrectPoly(Polymorphism correctPoly)
        {
            this.correctWeightPolys += correctPoly.getMutationRate();
            this.correctPolys.Add(correctPoly);
            this.unusedPolys.Remove(correctPoly);
        }
        public  void addExpectedPoly(Polymorphism correctPoly)
        {
            this.expectedWeightPolys += correctPoly.getMutationRate();
            this.allCheckedPolys.Add(correctPoly);
        }
        public  void removeExpectedPoly(Polymorphism currentPoly)
        {
            Polymorphism found = null;
            foreach (Polymorphism poly in this.allCheckedPolys)
            {
                if ((poly.Position == currentPoly.Position) && (poly.Mutation == currentPoly.Mutation))
                {
                    this.expectedWeightPolys -= this.allCheckedPolys[this.allCheckedPolys.IndexOf(poly)].getMutationRate();
                    found = poly;
                    Polymorphism newPoly = new Polymorphism(currentPoly);
                    newPoly.BackMutation = false;
                    this.correctedBackmutations.Add(new Polymorphism(newPoly));
                }
            }
            this.allCheckedPolys.Remove(found);
        }
        public  void removeCorrectPoly(Polymorphism currentPoly)
        {
            Polymorphism found = null;
            foreach (Polymorphism poly in this.correctPolys)
            {
                if ((poly.Position == currentPoly.Position) && (poly.Mutation == currentPoly.Mutation))
                {
                    this.correctWeightPolys -= correctPolys[this.correctPolys.IndexOf(poly)].getMutationRate();
                    if (!currentPoly.BackMutation)
                    {
                        this.unusedPolys.Add(currentPoly);
                    }
                    found = poly;
                    Polymorphism newPoly = new Polymorphism(currentPoly);
                    newPoly.BackMutation = false;
                    this.correctedBackmutations.Add(newPoly);
                }
            } this.correctPolys.Remove(found);
        }

        public  void addMissingOutOfRangePoly(Polymorphism correctPoly)
        {
            this.missingWeightPolys += correctPoly.getMutationRate();
            this.missingPolysOutOfRange.Add(correctPoly);
        }
        public  void removeMissingOutOfRangePoly(Polymorphism correctPoly)
        {
            this.missingWeightPolys -= correctPoly.getMutationRate();
            this.missingPolysOutOfRange.Add(correctPoly);
        }
        public  List<Polymorphism> UnusedNotInRange
        {
            set {
                this.unusedPolysNotInRange = value;
            }
        }
        public  int CompareTo(SearchResult o)
        {
            if (Rank > o.Rank)
            {
                return -1;
            }
            if (Rank < o.Rank)
            {
                return 1;
            } return 0;
        }
        public  PhyloTreePath UsedPath
        {
            get
            {
                return this.usedPath;
            }
        }
        public  void extendPath(PhyloTreeNode newNode)
        {
            this.usedPath.add(newNode);
        }
        public  List<Polymorphism> CorrectedBackmutations
        {
            get {
                return this.correctedBackmutations;
            }
        }
        public  PhyloTreePath PhyloTreePath
        {
            get {
                return this.usedPath;
            }
        }
    }
}