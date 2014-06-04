using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;
using Bio.Algorithms.Alignment;
using Bio.Algorithms.Assembly.Graph;
using Smrf.NodeXL.Adapters;
using Smrf.NodeXL.Core;
using MitoDataAssembler;
using Bio;

namespace MitoDataAssembler.GraphVisualization
{
    public class GraphGenerator
    {
        string outFileName;
        DeBruijnGraph _graph;
        public List<MegaNode> nodes = new List<MegaNode>();
        public GraphGenerator(string outPutFileName,DeBruijnGraph graph)
        {
            this.outFileName = outPutFileName;
            this._graph = graph;
        }
        /// <summary>
        /// Condense redundant paths down to simple paths
        /// </summary>
        /// <returns>List of simple paths.</returns>
        private IList<ISequence> CreateMegaNodes()
        {
            foreach(DeBruijnNode node in _graph.GetNodes())
            {
            IList<ISequence> paths = new List<ISequence>();
            Parallel.ForEach(this._graph.GetNodes(), node =>
            {
                int validLeftExtensionsCount = node.LeftExtensionNodesCount;
                int validRightExtensionsCount = node.RightExtensionNodesCount;

                if (validLeftExtensionsCount + validRightExtensionsCount == 0)
                {
                    // Island. Check coverage
                    if (Double.IsNaN(_coverageThreshold))
                    {
                        if (createContigSequences)
                        {
                            lock (paths)
                            {
                                paths.Add(_graph.GetNodeSequence(node));
                            }
                        }
                    }
                    else
                    {
                        if (node.KmerCount < _coverageThreshold)
                        {
                            node.MarkNodeForDelete();
                        }
                    }
                }
                else if (validLeftExtensionsCount == 1 && validRightExtensionsCount == 0)
                {
                    TraceSimplePath(paths, node, false, createContigSequences);
                }
                else if (validRightExtensionsCount == 1 && validLeftExtensionsCount == 0)
                {
                    TraceSimplePath(paths, node, true, createContigSequences);
                }
            });

            return paths;
        }

    }
}
