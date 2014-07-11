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
    /// Removes pathological sequences defined as:
    /// * All of one basepair type ("AAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTT", etc.)
    /// * Connected to themselves
    /// </summary>
    class PathologicalSequencePurger
    {        
	    /// <summary>
	    /// Delete nodes marked for erosion. Update adjacent nodes to update their extension tables.
	    /// </summary>
	    /// <param name="graph">De Bruijn Graph.</param>
	    public static int RemovePathologicalNodes(DeBruijnGraph graph)
	    {
		    //Basic strategy here, start at all reference nodes, go find everything that isn't in there
		    //and remove it.
		    DeBruijnGraph.ValidateGraph (graph);

            var badSeq = Enumerable.Repeat((byte)'A',graph.KmerLength).ToArray();
            var seq = new Bio.Sequence(Bio.Alphabets.DNA, badSeq,false);
            var badkmer1 = KmerData32.GetKmers(seq,graph.KmerLength).First().KmerData;

            badSeq = Enumerable.Repeat((byte)'G',graph.KmerLength).ToArray();
            seq = new Bio.Sequence(Bio.Alphabets.DNA, badSeq,false);
            var badkmer2 = KmerData32.GetKmers(seq,graph.KmerLength).First().KmerData;
            var badNodeCount = 0;
            foreach (var x in graph.GetNodes ()) {
				if (x.NodeValue.KmerData == badkmer1 ||
                    x.NodeValue.KmerData == badkmer2 ||
                    x.ContainsSelfReference) 
                {
					x.MarkNodeForDelete ();
                    Interlocked.Increment(ref badNodeCount);
				}
             }
			
            foreach(var node in graph.GetNodes())
            {
				node.RemoveMarkedExtensions();
			}
       
			//Now to delete them, since they are not connected to anything we are keeping,
			//no need to alter the graph structure
			graph.RemoveMarkedNodes ();
            return badNodeCount;
        }
	}
}