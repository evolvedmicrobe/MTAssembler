using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly;
using Bio.Algorithms.Assembly.Graph;
namespace MitoDataAssembler.Utils
{
    public static class CircularGenomeCaseHandlers
    {
        public readonly static int MT_Genome_Length = (int) HaploGrepSharp.ReferenceGenome.ReferenceSequence.Count;
        /// <summary>
        /// Decide if it makes more sense for this segment to "span" the origin
        /// or it makes more sense to think of it as going the other way around.
        /// </summary>
        /// <param name="putativeLeftSide">Putative left side</param>
        /// <param name="putatitiveRightSide">Putative right side</param>
        /// <param name="segmentLength">Segment length</param>
        /// <returns></returns>
        public static bool SegmentLoopsAround(double putativeLeftSide, double putatitiveRightSide, int segmentLength)
        {
            //space across span
            double right = Math.Min(putativeLeftSide, putatitiveRightSide);
            double left = Math.Max(putativeLeftSide, putatitiveRightSide);
            double side1 = Math.Abs(MT_Genome_Length - left) + right;
            //space on other side
            double side2 = Math.Abs(putatitiveRightSide - putativeLeftSide);
            side1 = Math.Abs(segmentLength - side1);
            side2 = Math.Abs(segmentLength - side2);
            if (side1 < side2)
                return true;
            else
                return false;
        }

        /// <summary>
        /// For a sequence ordered clockwise around the mtDNA, calculate the amount it spans on the reference, if loops, we adjust accordingly.
        /// </summary>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <param name="loops"></param>
        /// <returns></returns>
        public static int CalculateReferenceSpan(int start, int end, bool loops)
        {
            if (!loops)
            {
                return end - start + 1;
            }
            else
            {
                return (MT_Genome_Length - start +1 ) + end;
            }
        }

        /// <summary>
        /// Counts the number of times a node doesn't follow the node before or after it.
        /// </summary>
        /// <param name="nodes"></param>
        /// <returns></returns>
        public static int CountNumberOfTimesDirectionChanges(IEnumerable<DeBruijnNode> nodes)
        {
            var v=nodes.Where(x=>x.IsInReference).ToList();
            /* Now get a count of times the first node is after
             * the second node
             */
            if (v.Count < 2) 
            {
                return 0;
            }
            else
            {
                var difs = Enumerable.Zip(v.Skip(1),v.Take(v.Count-1),
                            (x,y)=> { 
                                if (x.ReferenceGenomePosition > y.ReferenceGenomePosition) 
                                    return 1;
                                else 
                                    return 0;}).Sum();
                int maxDifs = v.Count - 1;
                /* Depending on the most typical direction, the smallest 
                 * number of switches is the minimum of the following 
                 */
                return Math.Min(difs, maxDifs - difs);
            }
           
        }
    }
}
