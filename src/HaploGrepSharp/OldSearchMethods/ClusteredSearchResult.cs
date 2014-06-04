using System;
using System.Collections;
using System.Collections.Generic;
namespace HaploGrepSharp
{
    public class ClusteredSearchResult : IComparable<ClusteredSearchResult>
    {
        private int rankedPosition = -1;
        private List<SearchResult> cluster = new List<SearchResult>();
        public ClusteredSearchResult(int position)
        {
            this.rankedPosition = position;
        }
        public static List<ClusteredSearchResult> createClusteredSearchResult(List<SearchResult> unclusteredResults, Haplogroup hg)
        {
            List<ClusteredSearchResult> clusteredSearchResult = new List<ClusteredSearchResult>();
            //sort the results by top hit
            unclusteredResults.Sort();
            int i = 0;
            double currentRank = -1.0D;
            bool foundPredefinedHG = false;
            foreach (SearchResult currentResult in unclusteredResults)
            {
                if (i == 50){
                    break;
                } 
                if (currentRank != currentResult.Rank)
                {
                    clusteredSearchResult.Add(new ClusteredSearchResult(i + 1));
                    i++;
                    currentRank = currentResult.Rank;
                } 
                clusteredSearchResult[i - 1].Cluster.Add(currentResult);
                if (currentResult.Haplogroup.Equals(hg))
                {
                    foundPredefinedHG = true;
                }
            }
            if ((!foundPredefinedHG) && (!hg.ToString().Equals("")))
            {
                int pos = 1;
                foreach (SearchResult currentResult in unclusteredResults)
                {
                    if (currentResult.Haplogroup.Equals(hg))
                    {
                        clusteredSearchResult.Add(new ClusteredSearchResult(pos));
                        clusteredSearchResult[i].Cluster.Add(currentResult);
                        break;
                    } pos++;
                }
            }
            return clusteredSearchResult;
        }
        public int RankedPosition
        {
            get
            {
                return this.rankedPosition;
            }
        }
        public string Haplogroup
        {
            get
            {
                return Cluster[0].Haplogroup.ToString();
            }
        }
        public bool containsSuperhaplogroup(Haplogroup haplogroup)
        {
            if (haplogroup.isSuperHaplogroup(Cluster[0].Haplogroup))
            {
                return true;
            } return false;
        }
        public int CompareTo(ClusteredSearchResult o)
        {
            if (((SearchResult)Cluster[0]).Rank > ((SearchResult)o.Cluster[0]).Rank)
            {
                return -1;
            }
            if (((SearchResult)Cluster[0]).Rank < ((SearchResult)o.Cluster[0]).Rank)
            {
                return 1;
            } return 0;
        }
        public string ToString()
        {
            string result = "";
            System.Globalization.NumberFormatInfo df = new System.Globalization.NumberFormatInfo();//"0.000");
            result = result + "\t-------------------------------------------------------------------\n";
            foreach (SearchResult currentResult in Cluster)
            {
                result = result + "\tHaplogroup: " + currentResult.Haplogroup + "\n";
                result = result + "\tFinalRank: " + currentResult.Rank + " \n";
                result = result + "\tCorrect polys in test sample ratio: " + currentResult.CorrectPolyInTestSampleRatio + " (" + currentResult.CorrectWeightPolys.ToString() + " / " + currentResult.UsedWeightPolys.ToString() + ")" + "\n";
                result = result + "\tCorrect polys in haplogroup ratio: " + currentResult.CorrectPolyInHaplogroupRatio.ToString() + " (" + currentResult.CorrectWeightPolys.ToString() + " / " + currentResult.ExpectedWeightPolys.ToString() + ")" + "\n";
                result = result + "\t\tExpected\tCorrect\t\tUsed polys\tWeight\n";
                currentResult.CheckedPolys.Sort();
                ArrayList unusedPolys = new ArrayList();
                unusedPolys.AddRange(currentResult.Sample.Polymorphismn);
                foreach (Polymorphism current in currentResult.CheckedPolys)
                {
                    string fluctString = current.getMutationRate().ToString();
                    result = result + "\t\t" + current.ToString();
                    result = result + "\t\t";
                    if (currentResult.CorrectPolys.Contains(current))
                    {
                        result = result + current.ToString();
                    } result = result + "\t\t";
                    if (currentResult.Sample.Polymorphismn.Contains(current))
                    {
                        result = result + current.ToString();
                        unusedPolys.Remove(current);
                    }
                    result = result + "\t\t" + fluctString;
                    result = result + "\n";
                }
                foreach (Polymorphism current in unusedPolys)
                {
                    string fluctString = current.getMutationRate().ToString();
                    result = result + "\t\t\t\t\t\t" + current;
                    result = result + "\t\t" + fluctString + "\n";
                } result = result + "\t\t-------------------------------------------------------------\n";
                result = result + "\t\t" + currentResult.ExpectedWeightPolys.ToString() + "\t\t" + currentResult.CorrectWeightPolys.ToString() + "\t\t" + currentResult.UsedWeightPolys.ToString() + "\n";

                result = result + "\n\n";
            }
            return result;
        }
        public virtual List<SearchResult> Cluster
        {
            set
            {
                this.cluster = value;
            }
            get
            {
                return this.cluster;
            }
        }
        public virtual PhyloTreePath getPhyloTreePath(Haplogroup haplogroup)
        {
            foreach (SearchResult currentResult in this.cluster)
            {
                if (currentResult.Haplogroup.Equals(haplogroup))
                {
                    return currentResult.UsedPath;
                }
            }
            return null;
        }
        public virtual PhyloTreePath getPhyloTreePath(int index)
        {
            return ((SearchResult)this.cluster[index]).UsedPath;
        }
    }
}