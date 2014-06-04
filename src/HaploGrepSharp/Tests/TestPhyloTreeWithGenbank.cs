using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio;
using Bio.IO.FastA;
using HaploGrepSharp.TreeUtilities;

namespace HaploGrepSharp.Tests
{
    /// <summary>
    /// A class to load phylotree, load the downloaded accession numbers, and then verify that the mutations match.
    /// </summary>
    class TestPhyloTreeWithGenbankData
    {
        public static void RunTests()
        {
            Dictionary<string, Sequence> loadedGenomes = new Dictionary<string, Sequence>();
            System.IO.StreamWriter sw = new System.IO.StreamWriter("FailedTests.txt");
            sw.WriteLine("Haplogroup\tAccession\tReasonFailed\tInGenbankNotHaplogrep\tInHaplogrepNotGenbank");
            //First to load up the accession numbers
            string fileName = "PhyloTreeSequences.fasta.gz";
            FastAZippedParser fap = new FastAZippedParser(fileName);
            foreach (var seq in fap.Parse())
            {
                var accession = seq.ID.Split('|')[3];
                loadedGenomes[accession] = seq as Sequence;
                if(accession.Contains('.')) {
                    loadedGenomes[accession.Split('.')[0]]=seq as Sequence;
                }
            }
            //now to grab all the trees 
            var tree = TreeLoader.LoadTree(CONSTANTS.TREE_XML_FILE);
            int FailedTests = 0;
            int PassedTests = 0;
            int AmbiguousFail = 0;
            int SharedPolys = 0;
            int UnSharedPolys = 0;
            HashSet<int> BadPositions=new HashSet<int>();
            foreach (var node in tree.GetAllChildren())
            {
                //if (node.haplogroup.id != "J1c3g")
                //{ continue; }
                try
                {
                    //ignore empty and thousand genomes
                    if (!String.IsNullOrEmpty(node.haplogroup.AccessionId) &&
                        !node.haplogroup.AccessionId.StartsWith("NA")
                        && !node.haplogroup.AccessionId.StartsWith("HG"))
                    {
                        var id = node.haplogroup.AccessionId;
                        if (id == "NC_012920") { continue; }
                        
                        var seq = loadedGenomes[id];
                        if (seq.Alphabet.HasAmbiguity)
                        {
                            //sw.WriteLine(node.haplogroup.id + "\t" + node.haplogroup.AccessionId + "\tAmbiguousSequence");
                            AmbiguousFail++;
                            continue;
                        }
                        //now check for differences 
                        var origfoundPolys = NewSearchMethods.GenomeToHaploGrepConverter.FindPolymorphisms(seq).Where(x => !CONSTANTS.EXCLUDED_POSITIONS.Contains(x.position)).ToList();
                        var currentPolys = node.Mutations.ToList();//copy 
                        List<Polymorphism> removeLater = new List<Polymorphism>();
                        var foundPolys = origfoundPolys.ToList();
                        foreach (var poly in origfoundPolys)
                        {
                            if (currentPolys.Contains(poly))
                            { SharedPolys++; currentPolys.Remove(poly); foundPolys.Remove(poly); }
                        }
                        currentPolys.Sort((x, y) => x.position.CompareTo(y.position));
                        currentPolys = CONSTANTS.CommonPolymorphismFilter(currentPolys).ToList();
                        foundPolys = CONSTANTS.CommonPolymorphismFilter(foundPolys).ToList();
                        if (currentPolys.Count > 0)// || foundPolys.Count > 0)
                        {
                            UnSharedPolys += currentPolys.Count + foundPolys.Count;
                            sw.WriteLine(node.haplogroup.id
                                + "\t" + node.haplogroup.AccessionId
                                + "\tMismatchedPolymorphisms"
                                + "\t" + String.Join("-", foundPolys.Select(x => (object)x).ToArray())
                                + "\t" + String.Join("-", currentPolys.Select(x => (object)x).ToArray())
                                );
                            FailedTests++;
                            foreach (var v in foundPolys) { BadPositions.Add(v.position); }
                            foreach (var v in currentPolys) { BadPositions.Add(v.position); }
                        }
                        else { PassedTests++; }

                    }
                }
                catch (Exception thrown)
                {
                    sw.WriteLine(node.haplogroup.id
                                + "\t" + node.haplogroup.AccessionId
                                + "\t" + thrown.Message);                             
                }
            }
            sw.Close();
            sw = new System.IO.StreamWriter("TestReport.txt");
            sw.WriteLine("Failed Tests\t" + FailedTests.ToString());
            sw.WriteLine("Passed Tests\t" + PassedTests.ToString());
            sw.WriteLine("Ambiguous Sequences\t" + AmbiguousFail.ToString());
            sw.WriteLine("Polymorphisms overlapping between genbank and haplogrep\t" + SharedPolys.ToString());
            sw.WriteLine("Not Overlapping\t" + UnSharedPolys.ToString());
            sw.WriteLine("Total Discordant Positions\t" + BadPositions.Count.ToString());
            sw.WriteLine("Bad Position Listing");
            var bps = BadPositions.ToList();
            bps.Sort();
            foreach (var bp in bps)
            {
                sw.WriteLine(bp.ToString());
            }
            sw.Close();
        }

        public static void OutputHaplogroups()
        {
 
        }

    }
}
