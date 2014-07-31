using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly;
using Bio.Algorithms.Assembly.Padena;
using Bio.Algorithms.Assembly.Graph;
using System.Diagnostics;
using HaploGrepSharp;

namespace MitoDataAssembler
{
    using OutputColumn = OutputColumn<DeletionAnalysis>;
    /// <summary>
    /// Does not allow deletions to span the control region.
    /// </summary>
    public class DeletionAnalysis
    {
        
        private static int DeletionReportCounter=0;

        #region REPORT VALUES
        public int DeletionNumber{get;private set;}
        public bool ReferenceValuesChangeMonotonically
        {
            get;
            private set;
        }
        public bool HasDeletion;
        /// <summary>
        /// Does this assembly only have 2 sets of extensions with greater than one option?
        /// </summary>
        public bool SimpleBifurcation;

        /// <summary>
        /// If only two molecules are present in the entire assembly, than this gives the fractional split decided by looking at paths.
        /// </summary>
        public double AvgFractionBySplit;

        /// <summary>
        /// Size of gaps between adjacent nucmer alignments (separated by ";")
        /// </summary>
        public string SizeOfDeletionsSeen;

        public int StartReference;

        public int EndReference;

        public string DeletedRegions;
        #endregion
    
        public PossibleAssembly Assembly;

        public DeletionAnalysis(PossibleAssembly assemblyToCheck )
        {

            DeletionNumber=DeletionReportCounter++;
            Assembly = assemblyToCheck;
            LookForDeletion();
        }
        private void LookForDeletion()
        {
            bool movingUp;

            var v=Assembly.Where(x=>x.IsInReference).ToList();
            
            var difs=Enumerable.Zip(v.Skip(1),v.Take(v.Count-1),
                                    (x,y)=> { 
                                        if (x.ReferenceGenomePosition > y.ReferenceGenomePosition) 
                                            return 1;
                                        else 
                                            return 0;}).Sum();
            // If monotonically changing, should only not change once (when it goes around the circle).
            if (difs < 2 || difs > (v.Count - 3))
                ReferenceValuesChangeMonotonically = true;
            // Now which way is it increasing, up (big sum) or down (small sum)
            if (difs > (v.Count / 2) )
                movingUp = true;
            else
                movingUp = false;

            if (!movingUp)
            {
                Assembly.ReversePath();
                movingUp=true;
            }
            // Only report for sensible assemblies
            if (ReferenceValuesChangeMonotonically)
            {
                Assembly.FinalizeAndOrientToReference();
                // Get Alignments
                var alns = HaploGrepSharp.ReferenceGenome.GetDeltaAlignments(Assembly.Sequence).SelectMany(x => x).ToList();
                if (alns.Count > 0)
                {
                    StartReference = (int) alns.First().FirstSequenceStart;
                    EndReference = (int) alns.Last().FirstSequenceEnd;
                }
                SizeOfDeletionsSeen = String.Empty;
                if (alns.Count > 1)
                {
                    HasDeletion = true;
                    StringBuilder sb = new StringBuilder();
                    List<int> DeletionSizes = new List<int>();
                    for (int i = 0; i < (alns.Count - 1); i++)
                    {
                        var s = ReferenceGenome.ConvertTorCRSPosition((int)alns[i].FirstSequenceEnd);
                        var e = ReferenceGenome.ConvertTorCRSPosition((int)alns[i + 1].FirstSequenceStart);
                        sb.Append(s.ToString());
                        sb.Append("-");
                        sb.Append(e.ToString());
                        sb.Append(";");
                        DeletionSizes.Add(e - s + 1);
                    }
                    DeletedRegions = sb.ToString();
                    SizeOfDeletionsSeen = String.Join(";", DeletionSizes.Select(x => x.ToString()));
                }
                
                // Now see if we can get the fractional evidence for this.
                // Note this code can get very nasty as it we have multiple
                // nodes with 2, we can look for 
                var totBifurcations = 0;
               
                double avg = 0.0;
                var totL = Assembly.Where(x => x.LeftExtensionNodesCount == 2).ToList();
                var torR = Assembly.Where(x => x.RightExtensionNodesCount == (byte)2).ToList();
                for (int i = 0; i < (Assembly.Count -1); i++)
                {
                    var cnode = Assembly[i];
                    var lefts = cnode.GetLeftExtensionNodes().ToList();
                    var rights = cnode.GetRightExtensionNodes().ToList();
                    List<List<DeBruijnNode>> neighbors = new List<List<DeBruijnNode>>() { lefts, rights };
                    foreach (var neighbor in neighbors)
                    {
                        if (neighbor.Count == 2)
                        {
                            var tot = neighbor.Sum(z => (double)z.KmerCount);
                            var cur = (double)neighbor.Where(z => z == Assembly[i - 1] || z == Assembly[i + 1]).First().KmerCount;
                            avg += cur / tot;
                            totBifurcations++;
                        }
                        else if (neighbor.Count > 2)
                        {
                            totBifurcations = 100; //arbitrarily set to too high a value
                            break;
                        }
                    }
                }
                if (totBifurcations == 2)
                {
                    SimpleBifurcation = true;
                    avg *= .5;// .5 * (a + b) = Average

                }
                AvgFractionBySplit = SimpleBifurcation ? avg : Double.NaN;
                
            }            
        }
        public static string DeletionReportHeaderLine()
        {
            return String.Join(",", OutputColumnCollection.Select(x => x.Name).ToArray());
        }
        public IEnumerable<string> DeletionReportDataLines()
        {
            yield return String.Join(",", OutputColumnCollection.Select(x => x.GetValue(this)).ToArray());
        }

		/// <summary>
		/// Collection of columns to output with various summary statistics.
		/// </summary>
		public static List<OutputColumn> OutputColumnCollection = new List<OutputColumn>() {
            new OutputColumn("SequenceNumber", x => x.DeletionNumber.ToString()),
            new OutputColumn("HasDeletion", x => x.HasDeletion.ToString()),
            new OutputColumn("StartRef",x=>x.StartReference.ToString()),
            new OutputColumn("EndRef",x=>x.EndReference.ToString()),
            new OutputColumn("Deletion_Size(s)",x=>x.SizeOfDeletionsSeen),
            new OutputColumn("Deleted_Regions", x => x.DeletedRegions),
            new OutputColumn("SimpleBifurcation",x=>x.SimpleBifurcation.ToString()),
            new OutputColumn("FractionEvidence",x=>x.AvgFractionBySplit.ToString()),
            new OutputColumn("RefChangesMonotonically",x=>x.ReferenceValuesChangeMonotonically.ToString()),
            new OutputColumn("AvgCoverage",x=> x.Assembly.AvgKmerCoverage.ToString()),
            new OutputColumn("AssemblyLength",x=>x.Assembly.DirtySequence.Count.ToString()),
            new OutputColumn("SequenceOfAssembly",x=>x.Assembly.DirtySequence.ConvertToString())

		};
    }

   
}
