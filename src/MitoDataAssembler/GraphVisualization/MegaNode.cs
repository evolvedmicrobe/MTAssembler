using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly;
using Bio.Algorithms.Assembly.Graph;
using Bio;

namespace MitoDataAssembler.GraphVisualization
{
   
    /// <summary>
    /// A Node that represents a contiguous path in a graph
    /// (that is several subnodes) that occur with no divergence.
    /// </summary>
    public class MegaNode
    {
        List<DeBruijnNode> constituentNodes=new List<DeBruijnNode>();
        public ISequence currentSequence;
        public MegaNode()
        {
        }
    }
}
