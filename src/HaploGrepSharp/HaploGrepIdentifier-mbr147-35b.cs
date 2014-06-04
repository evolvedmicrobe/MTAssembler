using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Collections;

namespace HaploGrep
{
	/// <summary>
	/// A class which takes an HSD file and outputs possible haplotypes
	/// </summary>
    class HaploGrepIdentifier
    {
        Session session;
        List<string> lines = new List<string>();
        public HaploGrepIdentifier(string inFile,string outFile)
        {           
                string phylotree = "phylotree15.xml"; 
                string fluctrates = "fluctRates15.txt";
                session = new Session();
                FileInfo FI = new FileInfo(inFile);
                importData(FI);
                determineHG(phylotree, fluctrates);
                exportResults(outFile);
        }
        private void importData(FileInfo file)
        {
            StreamReader sr = new StreamReader(file.OpenRead());
            
            string line = "";
            sr.ReadLine();
            while ((line = sr.ReadLine()) != null)
            {
                lines.Add(line);
            }
            SampleFile newSampleFile = new SampleFile(lines);
            session.CurrentSampleFile = newSampleFile;
        }
        private void determineHG(string phylotree, string fluctrates)
        {
            HaploSearchManager h1 = new HaploSearchManager(phylotree, fluctrates); //PhylotreeInstance.Instance.getPhylotree(phylotree, fluctrates);
            HaploSearch haploSearch = new HaploSearch(h1);
            foreach (TestSample currentSample in session.CurrentSampleFile.TestSamples)
            {
                IList<ClusteredSearchResult> result = null;
                result = haploSearch.search(currentSample);
                haploSearch.addRecommendedHaplogroups(result, currentSample);
                session.setClassificationResults(currentSample.SampleID, result);
               
            }
        }
        private void exportResults(string outFilename)
        {
            StringBuilder result = new StringBuilder();
            var sampleCollection = session.CurrentSampleFile.TestSamples;
            sampleCollection.Sort();
            result.Append("SampleID\tRange\tHaplogroup\tQuality\tPolymorphisms\n");
            if (sampleCollection != null)
            {
                foreach (TestSample sample in sampleCollection)
                {
                    result.Append(sample.SampleID + "\t");
                    SampleRange range = sample.SampleRanges;
                    var startRange = range.Starts;
                    var endRange = range.Ends;
                    string resultRange = "";
                    for (int i = 0; i < startRange.Count; i++)
                    {
                        if (startRange[i]==endRange[i])
                        {
                            resultRange = resultRange + startRange[i] + ";";
                        }
                        else
                        {
                            resultRange = resultRange + startRange[i] + "-" + endRange[i] + ";";
                        }
                    } 
                    result.Append(resultRange + "\t");
                    result.Append(sample.RecognizedHaplogroup + "\t");
                    result.Append(sample.ResultQuality + "\t");
                    getPolysUsed(result, session.getClassificationResults(sample.SampleID));
                    getPolysNotUsed(result, session.getClassificationResults(sample.SampleID));
                    getPolysUnused(result, session.getClassificationResults(sample.SampleID));
                    result.Append("\n");
                }
            } 
            StreamWriter fileWriter = new StreamWriter(outFilename);
            fileWriter.Write(result.ToString());
            fileWriter.Close();
        }
        public static void getPolysNotUsed(StringBuilder buffer, IList<ClusteredSearchResult> a)
        {
            var correctPolys = ((SearchResult)((ClusteredSearchResult)a[0]).Cluster[0]).CorrectPolys;
            var checkedPolys = ((SearchResult)((ClusteredSearchResult)a[0]).Cluster[0]).CheckedPolys;
            correctPolys.Sort();
            checkedPolys.Sort();
            foreach (Polymorphism checkedPoly in checkedPolys)
            {
                if (!correctPolys.Contains(checkedPoly))
                {
                    buffer.Append(checkedPoly + " (no)\t");
                }
            }
        }
        public static void getPolysUsed(StringBuilder buffer, IList<ClusteredSearchResult> a)
        {
            var correctPolys = ((SearchResult)((ClusteredSearchResult)a[0]).Cluster[0]).CorrectPolys;
            var checkedPolys = ((SearchResult)((ClusteredSearchResult)a[0]).Cluster[0]).CheckedPolys;
            correctPolys.Sort();
            checkedPolys.Sort();
            foreach (Polymorphism checkedPoly in checkedPolys)
            {
                if (correctPolys.Contains(checkedPoly))
                {
                    buffer.Append(checkedPoly + " (yes)\t");
                }
            }
        }
        public static void getPolysUnused(StringBuilder buffer, IList<ClusteredSearchResult> a)
        {
            var unusedPolys = ((SearchResult)((ClusteredSearchResult)a[0]).Cluster[0]).UnusedPolys;
            var unusedPolysNotInRange = ((SearchResult)((ClusteredSearchResult)a[0]).Cluster[0]).UnusedPolysNotInRange;
            unusedPolys.Sort();
            foreach (Polymorphism unused in unusedPolys)
            {
                if (unused.getMutationRate(((SearchResult)((ClusteredSearchResult)a[0]).Cluster[0]).PhyloString) == 0.0D)
                {
                    if (unused.MTHotspot)
                    {
                        buffer.Append(unused + " (hotspot)\t");
                    }
                    else
                    {
                        buffer.Append(unused + " (globalPrivateMut)\t");
                    }
                }
                else if (unusedPolysNotInRange.Contains(unused))
                {
                    buffer.Append(unused + " (outOfRange)\t");
                }
                else
                {
                    buffer.Append(unused + " (localPrivateMut)\t");
                }
            }
        }	
    }
}
