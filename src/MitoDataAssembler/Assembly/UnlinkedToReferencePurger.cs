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
	/// Detects nodes that are in Cliques where a kmer from the reference sequence is not present, and then removes them.
	/// </summary>
	public class UnlinkedToReferencePurger 
	{
		string Name
				{
					get { return "Unlinked to reference purger"; }
				}
		string Description {
		get { return "Removes nodes thar are not connected to a kmer from the reference genome along any path."; }
	}

		/// <summary>
				/// Delete nodes marked for erosion. Update adjacent nodes to update their extension tables.
				/// </summary>
				/// <param name="graph">De Bruijn Graph.</param>
	    public void RemoveUnconnectedNodes(DeBruijnGraph graph, IEnumerable<DeBruijnNode> referenceNodes)
	{
		//Basic strategy here, start at all reference nodes, go find everything that isn't in there
		//and remove it.
		DeBruijnGraph.ValidateGraph (graph);
		//Mark all nodes as not visited
		graph.SetNodeVisitState (false);
		//Now visit everyone that is connected to the reference somehow
			//This loop should spend basically all its time on the first node
		foreach (DeBruijnNode node in referenceNodes) {
			if (node.IsVisited)
				continue;
			else {
				visitAllConnectedNodes (node);
			}
		}		
	    //Now mark any unvisited node for deletion.
		Parallel.ForEach (graph.GetNodes(), new ParallelOptions () { MaxDegreeOfParallelism=Environment.ProcessorCount }, x => {
			if (!x.IsVisited) {
				x.MarkNodeForDelete ();
			}});
        Parallel.ForEach(
               graph.GetNodes(),
               (node) =>
               {
                   node.RemoveMarkedExtensions();
               });
			//Now to delete them, since they are not connected to anything we are keeping,
			//no need to alter the graph structure
			graph.RemoveMarkedNodes ();
		}
		/// <summary>
		/// Visits all connected nodes, not caring about orientation here.  Hope the graph
		/// is not so big this leads to an OOM exception.
		/// </summary>
		/// <param name="startNode">Start node.</param>
    	private void visitAllConnectedNodes( DeBruijnNode startNode)
	{
			Stack<DeBruijnNode> toProcess = new Stack<DeBruijnNode> (16000);
			toProcess.Push(startNode);
			//Visit all nodes, avoid function all recursion with stack.
			do {
				DeBruijnNode next =toProcess.Pop();
				next.IsVisited=true;
				foreach (DeBruijnNode neighbor in next.GetExtensionNodes()) {
					if (neighbor.IsVisited)
						continue;
					else {
						toProcess.Push(neighbor);
					}
				}
			}while(toProcess.Count>0);		

		}
	}
}
