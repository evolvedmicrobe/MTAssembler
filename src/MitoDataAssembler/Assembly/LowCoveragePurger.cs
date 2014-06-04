using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly.Graph;
using Bio.Algorithms.Assembly.Padena;
using System.Collections.Concurrent;
using System.Threading;
using Bio.Algorithms.Kmer;
using Bio.Util;
using System.Diagnostics;
using System.Globalization;

namespace MitoDataAssembler
{
#if FALSE
    class LowCoverageAndUnPaintedPurger :IGraphErrorPurger
    {
       
        string Name
        {
            get { return "Low Coverage Unpainted Purger"; }
        }

        string Description
        {
            get { return "Removes graphical nodes below a given threshold unless they are painted as being a k-mer found in the reference genome"; }
        }

        IEnumerable<int> ErodeGraphEnds(DeBruijnGraph graph, int erosionThreshold)
        {

            DeBruijnGraph.ValidateGraph(graph);
            IEnumerable<int> toReturn = new SortedSet<int>();
            int TotalRemoved = 0;
            Parallel.ForEach(graph.GetNodes(), node =>
                {
                    if (node.KmerCount < erosionThreshold && node.ReferenceGenomePosition != DeBruijnNode.NO_GENOME_POSITION_SET_FLAG)
                    {
                        TotalRemoved++;
                        node.MarkNodeForDelete();
                    }
                });
                // Remove eroded nodes. In the out parameter, get the list of new 
                // end-points that was created by removing eroded nodes.
                RemoveErodedNodes(graph);
                return new int[]{TotalRemoved};
        }


        /// <summary>
        /// Delete nodes marked for erosion. Update adjacent nodes to update their extension tables.
        /// </summary>
        /// <param name="graph">De Bruijn Graph.</param>
        private void RemoveErodedNodes(DeBruijnGraph graph)
        {
            bool eroded = false;
            Parallel.ForEach(
                graph.GetNodes(),
                (node) =>
                {
                    if (node.IsMarkedForDelete)
                    {
                        node.IsDeleted = true;
                        eroded = true;
                    }
                });

            if (eroded)
            {
                Parallel.ForEach(
                graph.GetNodes(),
                (node) =>
                {
                    node.RemoveMarkedExtensions();
                });
            }
        }

       
		#region IGraphErrorPurger implementation

		public void RemoveErroneousNodes (DeBruijnGraph deBruijnGraph, DeBruijnPathList nodesList)
		{
			throw new NotImplementedException ();
		}

		public DeBruijnPathList DetectErroneousNodes (DeBruijnGraph deBruijnGraph)
		{

		}

		public int LengthThreshold {
			get {
				throw new NotImplementedException ();
			}
			set {
				throw new NotImplementedException ();
			}
		}

		#endregion

       
    }
#endif
}

