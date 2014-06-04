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
	/// <summary>
	/// Deletes any Kmer that only shows up once.
	/// </summary>
	public class ThresholdCoverageNodeRemover
	{
        string Name
        {
            get { return "Unlinked to reference purger"; }
        }
        string Description
        {
            get { return "Removes nodes thar are not connected to a kmer from the reference genome along any path."; }
        }
        /// <summary>
        /// Remove any node with coverage below this value
        /// </summary>
        public int CoverageCutOff = -1;

        public ThresholdCoverageNodeRemover(int cutoff)
        {
            this.CoverageCutOff = cutoff;
        }

	/// <summary>
	/// Delete nodes marked for erosion. Update adjacent nodes to update their extension tables.
	/// </summary>
	/// <param name="graph">De Bruijn Graph.</param>
	public void RemoveLowCoverageNodes(DeBruijnGraph graph)
	{
		//Basic strategy here, start at all reference nodes, go find everything that isn't in there
		//and remove it.
		DeBruijnGraph.ValidateGraph (graph);
		//Mark all nodes as not visited
		//Now visit everyone that is connected to the reference somehow
	    //Now mark any unvisited node for deletion.
			if (Bio.CrossPlatform.Environment.GetRunningPlatform () != Bio.CrossPlatform.Environment.Platform.Mac) {
				Parallel.ForEach (graph.GetNodes (), new ParallelOptions () { MaxDegreeOfParallelism=Environment.ProcessorCount }, x => {
					if (x.KmerCount < CoverageCutOff) {
						x.MarkNodeForDelete ();
					}});
				Parallel.ForEach(
					graph.GetNodes(),
					(node) =>
					{
					node.RemoveMarkedExtensions();
				});
			} else {
				foreach (var x in graph.GetNodes ()) {
					if (x.KmerCount < CoverageCutOff) {
						x.MarkNodeForDelete ();
					}
				}
			foreach(var node in
				graph.GetNodes()){
				node.RemoveMarkedExtensions();
			}  
		}
       
			//Now to delete them, since they are not connected to anything we are keeping,
			//no need to alter the graph structure
			graph.RemoveMarkedNodes ();
		}
	

	}
}
