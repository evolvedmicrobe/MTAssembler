﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly.Graph;
using Bio.Algorithms.Assembly.Padena;
using Bio.Algorithms.Alignment;
using Bio.Variant;
using Bio;
using Bio.SimilarityMatrices;

namespace MitoDataAssembler
{
	public class IndelPathCollection
	{
		public class IndelData
		{
            public double MinKmerCount;
			
            public Sequence Seq;
			
            public List<DeBruijnNode> NodesInPath;
			
            public int LikelyStart, LikelyEnd;

			public IndelData(DeBruijnPath originalPath, int kmerSize)
			{
				MinKmerCount = originalPath.PathNodes.Min(v => v.KmerCount);
                NodesInPath = originalPath.PathNodes.ToList();
				//Now to decide the rough start and end of this sequence, also need to orient it.
                SetStartAndEnd();
                var path = new DeBruijnPath(NodesInPath);
				Seq = path.ConvertToSequence(kmerSize);                
			}
            private void CheckStartEndForSanity()
            {
                if ((LikelyEnd - LikelyStart) > (200 + NodesInPath.Count) || (LikelyEnd < LikelyStart))
                {
                    throw new InvalidProgramException("Failed to accurately define indel region based on matches");
                }
            }
            /// <summary>
            /// This code sets the "region" the read aligns to based on the neighboring nodes
            /// in the assembly which have reference annotations.
            /// 
            /// The code is brutally long as it deals with separate cases in order to track having reference kmers on either side.
            /// </summary>
            private void SetStartAndEnd()
            {
                //First to get left and right, then flip if need be.
                var positions = Enumerable.Range(0, NodesInPath.Count).
                    Where(x => NodesInPath[x].IsInReference).
                    Select(z => new DistanceLocation() { Distance = z, RefGenomeLocation = NodesInPath[z].ReferenceGenomePosition})
                    .ToList();

                if (positions.Count > 1)
                {
                    //simplest option
                    var f = positions.First();
                    var e = positions.Last();
                    if (e.RefGenomeLocation < f.RefGenomeLocation)
                    {
                        NodesInPath.Reverse();
                        SetStartAndEnd();
                    }
                    else
                    {
                        LikelyStart = f.RefGenomeLocation - f.Distance;
                        LikelyEnd = e.RefGenomeLocation + (NodesInPath.Count - 1 - e.Distance);
                        CheckStartEndForSanity();
                    }
                }
                else 
                {
                    var leftStart = GetLeftSide();
                    var rightStart = GetRightSide();
                    if (positions.Count == 1)
                    {
                        var curVal = positions[0];
                        if (leftStart != null)
                        {
                            var alt = leftStart.RefGenomeLocation;
                            if (alt > curVal.RefGenomeLocation)
                            {
                                NodesInPath.Reverse();
                                SetStartAndEnd();
                            }
                            else
                            {
                                LikelyStart = leftStart.RefGenomeLocation + leftStart.Distance;
                                LikelyEnd = curVal.RefGenomeLocation + (NodesInPath.Count - 1 - curVal.Distance);
                                CheckStartEndForSanity();
                            }
                        }
                        else if (rightStart != null)
                        {
                            var alt = rightStart.RefGenomeLocation;
                            if (alt < curVal.RefGenomeLocation)
                            {
                                NodesInPath.Reverse();
                                SetStartAndEnd();
                            }
                            else
                            {
                                LikelyEnd = rightStart.RefGenomeLocation - rightStart.Distance;
                                LikelyStart = curVal.RefGenomeLocation - curVal.Distance;
                                CheckStartEndForSanity();
                            }
                        }
                    }
                    else
                    {
                        if (rightStart == null || leftStart == null)
                        {
                            throw new InvalidProgramException("Left and right Starts did not both exist, cannot place indel.");
                        }
                        if (rightStart.RefGenomeLocation < leftStart.RefGenomeLocation)
                        {
                            NodesInPath.Reverse();
                            SetStartAndEnd();
                        }
                        else
                        {
                            LikelyStart = leftStart.RefGenomeLocation + leftStart.Distance;
                            LikelyEnd = rightStart.RefGenomeLocation - rightStart.Distance;
                            CheckStartEndForSanity();
                        }
                    }
                }
            }

            /// <summary>
            /// Follows all paths leaving a node and returns the first one that has a reference genome location.
            /// </summary>
            /// <param name="node"></param>
            /// <param name="grabRightSide"></param>
            /// <param name="curDistance">How far removed from node we are at present</param>
            /// <returns></returns>
            public DistanceLocation FollowNode(DeBruijnNode node, bool grabRightSide, int curDistance)
            {
                 var nextNodes =
                    grabRightSide ? node.GetRightExtensionNodesWithOrientation() : node.GetLeftExtensionNodesWithOrientation();
                 foreach (var neighbor in nextNodes)
                 {
                     if (neighbor.Key.IsInReference)
                     {
                         return new DistanceLocation() { Distance = curDistance, RefGenomeLocation = neighbor.Key.ReferenceGenomePosition };
                     }
                     else
                     {
                         var nextSideRight = !(neighbor.Value ^ grabRightSide);
                         var res = FollowNode(neighbor.Key, nextSideRight, curDistance + 1);
                         if (res != null)
                         {
                             return res;
                         }
                     }
                 }
                 return null;
            }
           
            public DistanceLocation GetLeftSide()
            {            
                //  For the first node we add its sequence, and the neighbor could be on the left or right side
                var cur_node = NodesInPath[0];
                var leftNodes = cur_node.GetLeftExtensionNodesWithOrientation().Where(x => x.Key == NodesInPath[1]).ToList();
                var rightNodes = cur_node.GetRightExtensionNodesWithOrientation().Where(x => x.Key == NodesInPath[1]).ToList();
                var goingLeft = leftNodes.Count == 1;
                if (!((leftNodes.Count == 1) ^ (rightNodes.Count == 1)))
                {
                    throw new InvalidProgramException();
                }
                return FollowNode(cur_node, goingLeft, 0);
            }

            public DistanceLocation GetRightSide()
            {
                //  For the first node we add its sequence, and the neighbor could be on the left or right side
                var last_node = NodesInPath.Last();
                var penultimate_node = NodesInPath[NodesInPath.Count - 2];
                var leftNodes = penultimate_node.GetLeftExtensionNodesWithOrientation().Where(x => x.Key == last_node).ToList();
                var rightNodes = penultimate_node.GetRightExtensionNodesWithOrientation().Where(x => x.Key == last_node).ToList();
                var goingToLeft = leftNodes.Count == 1;
                if (!((leftNodes.Count == 1) ^ (rightNodes.Count == 1)))
                {
                    throw new InvalidProgramException();
                }
                var next_direction = goingToLeft ? leftNodes.First() : rightNodes.First();
                var nextSideRight = !(next_direction.Value ^ !goingToLeft);
                return FollowNode(last_node, nextSideRight, 0);
            }



		}
        /// <summary>
        /// A struct that is useful for assessing where a node is.
        /// </summary>
        public class DistanceLocation
        {
            /// <summary>
            /// How far is this node from the start node?
            /// </summary>
            public int Distance;

            /// <summary>
            /// What is the RefGenomeLocation at that node?
            /// </summary>
            public int RefGenomeLocation;
        }

        public class IndelLocation : IComparable<IndelLocation>, IEquatable<IndelLocation>
        {
            public bool DeletionOnReference;
            public string InsertedSequence;
            public int Start, Length;
            public IndelLocation(bool deletionOnReference, int start, int length)
            {
                DeletionOnReference = deletionOnReference;
                Start = start;
                Length = length;
            }

            public override bool Equals(object obj)
            {
                var no = obj as IndelLocation;
                if (no != null)
                {
                    return this.Equals(no);
                }
                return false;
            }
            public override int GetHashCode()
            {
                var str = String.Join("-", DeletionOnReference.ToString(), InsertedSequence, Start.ToString(), Length.ToString());
                return str.GetHashCode();
            }
            public bool Equals(IndelLocation other)
            {
                return DeletionOnReference == other.DeletionOnReference &&
                    Start == other.Start &&
                    Length == other.Length;
            }



            public int CompareTo(IndelLocation other)
            {
                if (other == this || this.Equals(other))
                {
                    return 0;
                }
                else
                {
                    if (other.Start < this.Start)
                        return 1;
                    else if (other.Start > this.Start)
                        return -1;
                    else
                        return this.Length.CompareTo(other.Length);
                }
            }

            bool IEquatable<IndelLocation>.Equals(IndelLocation other)
            {
                return Equals(other);
            }
        }

        public const int AlignmentPadding = 5;

		public static List<ContinuousFrequencyIndelGenotype> CallIndelsFromPathCollection (DeBruijnPathList paths, DeBruijnGraph graph)
		{
			var sequences = paths.Paths.Select(z => new IndelData(z, graph.KmerLength)).ToList();
            var regionStart = sequences.Min(x => x.LikelyStart) - IndelPathCollection.AlignmentPadding - graph.KmerLength;
            var regionEnd = sequences.Max(x => x.LikelyEnd) + IndelPathCollection.AlignmentPadding + graph.KmerLength;
            var reference = HaploGrepSharp.ReferenceGenome.GetReferenceSequenceSection(regionStart, regionEnd);

            var algo = new Bio.Algorithms.Alignment.SmithWatermanAligner();

            // Setup the aligner with appropriate parameters
            algo.SimilarityMatrix = new Bio.SimilarityMatrices.DiagonalSimilarityMatrix (1, -1);
            algo.GapOpenCost = -2;
            algo.GapExtensionCost = -1;

            // Execute the alignment.
            //Now to go through and generate variants.
            Dictionary<IndelData, List<IndelLocation>> indels = new Dictionary<IndelData, List<IndelLocation>>();
            System.IO.StreamWriter sw = new System.IO.StreamWriter("test.txt");
            foreach (var s in sequences)
            {
                //Note, do not change alignment order here.
                var aln = algo.Align(reference.Seq, s.Seq);
                var res = aln[0].PairwiseAlignedSequences[0];
                sw.WriteLine(reference.Seq.ConvertToString());
                sw.WriteLine(s.Seq.ConvertToString());
                sw.WriteLine(res.ToString());
                sw.WriteLine();
                Console.WriteLine(reference.Seq.ConvertToString());
                Console.WriteLine(s.Seq.ConvertToString());
                Console.WriteLine(res.ToString());
                var indels_locs = FindIndels(res);
                indels[s] = indels_locs;
            }

            sw.Close();
            //Now to get unique starts and run them.  
            var locations = indels.Values.SelectMany(z => z).ToList().Distinct().GroupBy(z => z.Start);
            var toReturn = new List<ContinuousFrequencyIndelGenotype>(10);
            //typcially there will only be one of these
            foreach (var g in locations)
            {
                var g2 = g.ToList();
                var lspots = g2.Select(x=>x.DeletionOnReference).Distinct().Count();
                if(lspots > 1)
                {
                    throw new NotImplementedException("Same location had indels present on both reference and reads.  This is an edge case not handled");
                }
                //add a fake null.
                var first = g2[0];
                var no_indel = new IndelLocation(first.DeletionOnReference, first.Start, 0);
                no_indel.InsertedSequence = String.Empty;
                g2.Add(no_indel);
                //sort by location
                g2.Sort();
                double[] counts = new double[g2.Count];
                //now add counts from each.
                foreach (var s in sequences)
                {
                    var cur = indels[s];
                    int index = 0;
                    var atLoc = cur.Where(x => x.Start == first.Start).FirstOrDefault();
                    if (atLoc == null)
                    {
                        counts[index]+=s.MinKmerCount;
                    }
                    else 
                    {
                        bool found = false;
                        for (int i = 1; i < g2.Count; i++)
                        {
                            if (g2[i].Equals(atLoc))
                            {
                                counts[i] += s.MinKmerCount;
                                found = true;
                                break;
                            }
                        }
                        if (!found)
                        {
                            throw new InvalidProgramException("Point should never be reached");
                        }
                    }
                }
                var types = g2.Select(x => x.InsertedSequence).ToList();
                var ref_loc = HaploGrepSharp.ReferenceGenome.ConvertTorCRSPosition(first.Start + regionStart);
                var indel_call = new ContinuousFrequencyIndelGenotype(first.DeletionOnReference, types, counts, ref_loc);
                toReturn.Add(indel_call);
            }
            return toReturn;
		}

        private static List<IndelLocation> LeftAlignIndels(PairwiseAlignedSequence aln, List<IndelLocation> curIndels)
        {
            var haplo = aln.SecondSequence;
            var refer = aln.FirstSequence;
            var toReturn = new List<IndelLocation>();
            foreach (var indel in curIndels)
            {
                int start = indel.Start;
                int len = indel.Length;
                ISequence deleted, notDeleted;
                if (indel.DeletionOnReference)
                {
                    deleted = refer;
                    notDeleted = haplo;
                }
                else {
                    deleted = haplo;
                    notDeleted = refer;
                }
                while ((start+len) < notDeleted.Count && start >0)
                {
                    //is a ,a
                    bool matchBefore = notDeleted[start - 1] == deleted [start - 1];
                    bool matchAfterShift = deleted[start - 1] == notDeleted[start + len - 1];
                    if (matchAfterShift && matchBefore)
                    {
                        indel.Start--;
                        start--;
                    }
                    else { break; }
                }
            }
            return curIndels;

        }

        private static List<IndelLocation> FindIndels(PairwiseAlignedSequence aln)
        {
            //Go through and find indels
            var haplo = aln.SecondSequence;
            var refer = aln.FirstSequence;
            var toReturn = new List<IndelLocation>();
            for (int i = 0; i < haplo.Count; i++)
            {
                if (haplo[i] == '-')
                {
                    int length = 1;
                    int start = i;
                    
                    do
                    {
                        if (haplo[++i] == '-')
                        {
                            length++;
                        }
                        else { break; }
                    } while (i < haplo.Count);
                    toReturn.Add(new IndelLocation(false, start, length));
                }
                if (refer[i] == '-')
                {
                    int length = 1;
                    int start = i;
                    do
                    {
                        if (refer[++i] == '-')
                        {
                            length++;
                        }
                        else { break; }
                    } while (i < haplo.Count);
                    
                    toReturn.Add(new IndelLocation(true, start, length));
                }
            }
            //left align them
            var lAlign = LeftAlignIndels(aln,toReturn);
            //now reposition them and get the missing sequence
            //be sure to account for any offsets introduced by past indels
            var refGapsSeen = 0;
            foreach(var laln in lAlign)
            {

                var relSeq = laln.DeletionOnReference ? haplo : refer;
                laln.InsertedSequence = new String(relSeq.Skip(laln.Start).Take(laln.Length).Select(x => (char)x).ToArray());
                laln.Start = laln.Start + (int)aln.SecondOffset - refGapsSeen;
                if (laln.DeletionOnReference) { refGapsSeen += laln.Length; }
            }
            return lAlign;
        }

	}
}
