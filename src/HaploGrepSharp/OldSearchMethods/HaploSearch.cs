using System.Collections;
using System.Collections.Generic;
using System.Xml;

namespace HaploGrepSharp
{
    public class HaploSearch
    {
        internal HaploSearchManager searchManager = null;
        public HaploSearch(HaploSearchManager searchManager) {
            this.searchManager = searchManager;
        }
        public IList<ClusteredSearchResult> search(TestSample testSample) {
            ///Create a list of all possible results, visiting each node.
            List<SearchResult> results = searchPhylotreeWrapper(testSample);
            //now cluster the results
            var clusteredResult = ClusteredSearchResult.createClusteredSearchResult(results, testSample.PredefiniedHaplogroup);
            results.Clear();
            return clusteredResult;
        }
        public void addRecommendedHaplogroups(IList<ClusteredSearchResult> result, TestSample sample) {
            sample.RecognizedHaplogroup = new Haplogroup(result[0].Haplogroup);
            double firstRank = result[0].Cluster[0].Rank * 100.0D;
            if (firstRank > 0.0D) {
                //decimal myDec = new decimal(firstRank);
                //myDec = myDec.setScale(1, 4);
                sample.ResultQuality = firstRank;//(double)myDec;
            }
            if (sample.PredefiniedHaplogroup.Equals(sample.RecognizedHaplogroup)) {
                sample.State = "identical";
            }
            else if ((sample.PredefiniedHaplogroup.isSuperHaplogroup(sample.RecognizedHaplogroup)) || (sample.RecognizedHaplogroup.isSuperHaplogroup(sample.PredefiniedHaplogroup))) {
                sample.State = "similar";
            }
            else {
                sample.State = "mismatch";
            }
        }
        private List<SearchResult> searchPhylotreeWrapper(TestSample sample) {
            var results = new List<SearchResult>();
            SearchResult rootResult = new SearchResult("rCRS, NC_012920", sample.Sample);
            XmlNode node = this.searchManager.PhyloTree.DocumentElement.SelectSingleNode("/phylotree");//.RootElement;
            searchPhylotree(node, results, sample, rootResult);
            return results;
        }
        /// <summary>
        /// Recursive function to search the XML tree and create a result at each node
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="results"></param>
        /// <param name="sample"></param>
        /// <param name="parentResult"></param>
        private void searchPhylotree(XmlNode parent, List<SearchResult> results, TestSample sample, SearchResult parentResult)
        {
            var children = parent.SelectNodes("haplogroup");
            foreach (XmlNode currentElement in children) {
                //string,string,parent result
                SearchResult newResult = new SearchResult(currentElement.Attributes.GetNamedItem("name").Value, parentResult);
                var details = currentElement.SelectSingleNode("details");
                var polys = details.SelectNodes("poly");
                PhyloTreeNode newNode = new PhyloTreeNode(new Haplogroup(currentElement.Attributes.GetNamedItem("name").Value));
                foreach (XmlNode currentPolyElement in polys)
                {
                    if (currentPolyElement.InnerText.Contains("X"))
                    { 
                        //System.Console.WriteLine("Skipping: " + currentPolyElement.InnerText); 
                        continue;
                    }
                    Polymorphism currentPoly = new Polymorphism(currentPolyElement.InnerText);
                    if (sample.SampleRanges.contains(currentPoly))
                    {
                        if (currentPoly.BackMutation) {
                            newResult.removeExpectedPoly(currentPoly);
                            newResult.removeCorrectPoly(currentPoly);
                            newNode.removeExpectedPoly(currentPoly);
                            newNode.removeCorrectPoly(currentPoly);
                            newNode.addExpectedPoly(currentPoly);
                            Polymorphism newPoly = new Polymorphism(currentPoly);
                            newPoly.BackMutation = false;
                            if (!newResult.Sample.contains(newPoly)) {
                                newNode.addCorrectPoly(currentPoly);
                            }
                        }
                        else if (newResult.Sample.contains(currentPoly)) {
                            newResult.addExpectedPoly(currentPoly);
                            newResult.addCorrectPoly(currentPoly);
                            newNode.addExpectedPoly(currentPoly);
                            newNode.addCorrectPoly(currentPoly);
                        }
                        else {
                            if (currentPoly.BackMutation) {
                                newResult.removeMissingOutOfRangePoly(currentPoly);
                            }
                            newResult.addExpectedPoly(currentPoly);
                            newNode.addExpectedPoly(currentPoly);
                        }
                    }
                    else {
                        newResult.addMissingOutOfRangePoly(currentPoly);
                        newNode.addNotInRangePoly(currentPoly);
                    }
                }
                newResult.UnusedNotInRange = sample.PolyNotinRange;
                results.Add(newResult);
                newResult.extendPath(newNode);
                searchPhylotree(currentElement, results, sample, newResult);
            }
        }
    }
}