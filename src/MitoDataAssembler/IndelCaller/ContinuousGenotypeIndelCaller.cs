using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly.Graph;
using Bio.Algorithms.Assembly.Padena;
namespace MitoDataAssembler.IndelCaller
{
    /// <summary>
    /// Calls indels by looking for redundant paths of differing lengths.
    /// </summary>
    public class ContinuousGenotypeIndelCaller 
    {

        #region Fields, Constructor, Properties

        /// <summary>
        /// Holds reference to assembler graph.
        /// </summary>
        private DeBruijnGraph graph;

        /// <summary>
        /// Initializes a new instance of the RedundantPathsPurger class.
        /// Takes user parameter for threshold. 
        /// </summary>
        /// <param name="length">Threshold length.</param>
        public ContinuousGenotypeIndelCaller (int length)
        {
            this.pathLengthThreshold = length;
        }

        /// <summary>
        /// Threshold for length of redundant paths.
        /// </summary>
        private int pathLengthThreshold;


        /// <summary>
        /// Gets or sets threshold for length of redundant paths.
        /// </summary>
        public int LengthThreshold
        {
            get { return this.pathLengthThreshold; }
            set { this.pathLengthThreshold = value; }
        }
        #endregion

        /// <summary>
        /// Detect nodes that are on redundant paths. 
        /// Start from any node that has ambiguous (more than one) extensions.
        /// From this node, trace path for each extension until either they 
        /// converge to a single node or threshold length is exceeded. 
        /// In case they converge, we have a set of redundant paths. 
        /// We pick the best path based on the kmer counts of the path nodes.
        /// All paths other than the best one are returned for removal.
        /// Locks: Method only does reads. No locking necessary here or its callees. 
        /// </summary>
        /// <param name="deBruijnGraph">De Bruijn Graph.</param>
        /// <returns>List of path nodes to be deleted.</returns>
        public DeBruijnPathList DetectIndels(DeBruijnGraph deBruijnGraph)
        {
            if (deBruijnGraph == null)
            {
                throw new ArgumentNullException("deBruijnGraph");
            }

            this.graph = deBruijnGraph;

            List<DeBruijnPathList> redundantPaths = new List<DeBruijnPathList>();
            Parallel.ForEach(
                deBruijnGraph.GetNodes(),
                node =>
                {
                    // Need to check for both left and right extensions for ambiguity.
                    if (node.RightExtensionNodesCount > 1)
                    {
                        TraceDivergingExtensionPaths(node, node.GetRightExtensionNodesWithOrientation(), true, redundantPaths);
                    }

                    if (node.LeftExtensionNodesCount > 1)
                    {
                        TraceDivergingExtensionPaths(node, node.GetLeftExtensionNodesWithOrientation(), false, redundantPaths);
                    }
                });

            var counts = redundantPaths.Select(x => String.Join(",", x.Paths.Select(z => z.PathNodes.Count).Distinct().Select(k=>k.ToString()))).ToList();
            var indelPaths = redundantPaths.Where(x => x.Paths.Select(z => z.PathNodes.Count).Distinct().Count()!=1).ToList();
            foreach (var v in indelPaths)
            {
                Console.WriteLine();
                Console.WriteLine((new DeBruijnPath(v.Paths.First().PathNodes.First())).ConvertToSequence(19).ConvertToString());
                Console.WriteLine((new DeBruijnPath(v.Paths.First().PathNodes.Last())).ConvertToSequence(19).ConvertToString());
                foreach (var k in v.Paths)
                {
                    Console.WriteLine(k.ConvertToSequence(19).ConvertToString());
                }
            }
            redundantPaths = RemoveDuplicates(redundantPaths);
            return DetachBestPath(redundantPaths);
        }

        /// <summary>
        /// Removes nodes that are part of redundant paths. 
        /// </summary>
        /// <param name="deBruijnGraph">De Bruijn graph.</param>
        /// <param name="nodesList">Path nodes to be deleted.</param>
        public void RemoveErroneousNodes(DeBruijnGraph deBruijnGraph, DeBruijnPathList nodesList)
        {
            if (this.graph == null)
            {
                throw new ArgumentNullException("deBruijnGraph");
            }

            DeBruijnGraph.ValidateGraph(deBruijnGraph);

            if (nodesList == null)
            {
                throw new ArgumentNullException("nodesList");
            }

            this.graph = deBruijnGraph;

            // Neighbors of all nodes have to be updated.
            HashSet<DeBruijnNode> deleteNodes = new HashSet<DeBruijnNode>(
                nodesList.Paths.AsParallel().SelectMany(nl => nl.PathNodes));

            // Update extensions for deletion
            // No need for read-write lock as deleteNode's dictionary is being read, 
            // and only other graph node's dictionaries are updated.
            Parallel.ForEach(
                deleteNodes,
                node =>
                {
                    foreach (DeBruijnNode extension in node.GetExtensionNodes())
                    {
                        // If the neighbor is also to be deleted, there is no use of updation in that case
                        if (!deleteNodes.Contains(extension))
                        {
                            extension.RemoveExtensionThreadSafe(node);
                        }
                    }
                });

            // Delete nodes from graph
            this.graph.RemoveNodes(deleteNodes);
        }

        /// <summary>
        /// Extract best path from the list of paths in each cluster.
        /// Take off the best path from list and return rest of the paths
        /// for removal.
        /// </summary>
        /// <param name="pathClusters">List of path clusters.</param>
        /// <returns>List of path nodes to be removed.</returns>
        private static DeBruijnPathList DetachBestPath(List<DeBruijnPathList> pathClusters)
        {
            return new DeBruijnPathList(
                pathClusters.AsParallel().SelectMany(paths => ExtractBestPath(paths).Paths));
        }

        /// <summary>
        /// Extract best path from list of paths. For the current cluster 
        /// of paths, return only those that should be removed.
        /// </summary>
        /// <param name="divergingPaths">List of redundant paths.</param>
        /// <returns>List of paths nodes to be deleted.</returns>
        private static DeBruijnPathList ExtractBestPath(DeBruijnPathList divergingPaths)
        {
            // Find "best" path. Except for best path, return rest for removal 
            int bestPathIndex = GetBestPath(divergingPaths);

            DeBruijnPath bestPath = divergingPaths.Paths[bestPathIndex];
            divergingPaths.Paths.RemoveAt(bestPathIndex);

            // There can be overlap between redundant paths.
            // Remove path nodes that occur in best path
            foreach (var path in divergingPaths.Paths)
            {
                path.RemoveAll(n => bestPath.PathNodes.Contains(n));
            }

            return divergingPaths;
        }

        /// <summary>
        /// Gets the best path from the list of diverging paths.
        /// Path that has maximum sum of 'count' of belonging k-mers is best.
        /// In case there are multiple 'best' paths, we arbitrarily return one of them.
        /// </summary>
        /// <param name="divergingPaths">List of diverging paths.</param>
        /// <returns>Index of the best path.</returns>
        private static int GetBestPath(DeBruijnPathList divergingPaths)
        {
            // We find the index of the 'best' path.
            long max = -1;
            int maxIndex = -1;

            // Path that has the maximum sum of 'count' of belonging k-mers is the winner
            for (int i = 0; i < divergingPaths.Paths.Count; i++)
            {
                long sum = divergingPaths.Paths[i].PathNodes.Sum(n => n.KmerCount);
                if (sum > max)
                {
                    max = sum;
                    maxIndex = i;
                }
            }

            return maxIndex;
        }

        /// <summary>
        /// Gets start node of redundant path cluster
        /// All paths in input are part of a redundant path cluster
        /// So all of them have the same start and the end node.
        /// Return the first node of first path.
        /// </summary>
        /// <param name="paths">List of redundant paths.</param>
        /// <returns>Start node of redundant path cluster.</returns>
        private static DeBruijnNode GetStartNode(DeBruijnPathList paths)
        {
            return paths.Paths.First().PathNodes.First();
        }

        /// <summary>
        /// Gets end node of redundant path cluster
        /// All paths in input are part of a redundant path cluster
        /// So all of them have the same start and the end node.
        /// Return the last node of first path.
        /// </summary>
        /// <param name="paths">List of redundant paths.</param>
        /// <returns>End node of redundant path cluster.</returns>
        private static DeBruijnNode GetEndNode(DeBruijnPathList paths)
        {
            return paths.Paths.First().PathNodes.Last();
        }

        /// <summary>
        /// Some set of paths will appear twice, one traced in forward direction
        /// and other in opposite. This method eliminate duplicates.
        /// </summary>
        /// <param name="redundantPathClusters">List of path cluster.</param>
        /// <returns>List of unique path clusters.</returns>
        private static List<DeBruijnPathList> RemoveDuplicates(List<DeBruijnPathList> redundantPathClusters)
        {
            // Divide the list into two groups. One with paths that do not 
            // have duplicates, and one with paths that do not have duplicate
            List<IGrouping<bool, DeBruijnPathList>> uniqueAndDuplicatedPaths =
            redundantPathClusters.AsParallel().GroupBy(pc1 =>
                redundantPathClusters.Any(pc2 =>
                    GetStartNode(pc1) == GetEndNode(pc2) && GetEndNode(pc1) == GetStartNode(pc2))).ToList();

            List<DeBruijnPathList> uniquePaths = new List<DeBruijnPathList>();
            foreach (IGrouping<bool, DeBruijnPathList> group in uniqueAndDuplicatedPaths)
            {
                if (!group.Key)
                {
                    // Add all paths that do have duplicates to final list
                    uniquePaths.AddRange(group);
                }
                else
                {
                    // Each element in this list contains a duplicate in the list
                    // Add only those where the start node has a sequence that is
                    // lexicographically greater than the end node sequence. This
                    // operation will eliminate duplicates effectively.
                    uniquePaths.AddRange(
                        group.AsParallel().Where(pc =>
                                GetStartNode(pc).NodeValue.CompareTo(
                                GetEndNode(pc).NodeValue) >= 0));
                }
            }

            return uniquePaths;
        }

        /// <summary>
        /// Traces diverging paths in given direction.
        /// For each path in the set of diverging paths, extend path by one node
        /// at a time. Continue this until all diverging paths converge to a 
        /// single node or length threshold is exceeded.
        /// If paths converge, add path cluster containing list of redundant 
        /// path nodes to list of redundant paths and return.
        /// </summary>
        /// <param name="startNode">Node at starting point of divergence.</param>
        /// <param name="divergingNodes">List of diverging nodes.</param>
        /// <param name="isRightExtension">Bool indicating direction of divergence.</param>
        /// <param name="redundantPaths">List of redundant paths.</param>
        private void TraceDivergingExtensionPaths(
            DeBruijnNode startNode,
            Dictionary<DeBruijnNode, bool> divergingNodes,
            bool isRightExtension,
            List<DeBruijnPathList> redundantPaths)
        {
            //maka a new path with each having the same start node, and a differing second node based on orientation
            List<PathWithOrientation> divergingPaths = new List<PathWithOrientation>(
                divergingNodes.Select(n =>
                    new PathWithOrientation(startNode, n.Key, n.Value)));
            int maxDivergingPathLength = 2;

            // Extend paths till length threshold is exceeded.
            // In case paths coverge within threshold, we break out of while.
            
            // Make a list of paths that we have finished following, this would be any path that
            // has reached a divergent end, possibly before the others have. 
            var finishedPaths = new List<PathWithOrientation>(divergingPaths.Count);
            while (maxDivergingPathLength <= this.pathLengthThreshold && finishedPaths.Count != divergingPaths.Count)
            {
                // Extend each path in cluster. While performing path extension 
                // also keep track of whether they have converged
                foreach (PathWithOrientation path in divergingPaths)
                {
                    if (finishedPaths.Contains(path))
                    {
                        continue;
                    }
                    /* We go left if we are already heading left in the same orientation, or if
                       we are heading right with a different orientation */
                    var grabLeftNext = isRightExtension ^ path.IsSameOrientation;
                    var endNode = path.Nodes.Last();
                    var nextNodes = grabLeftNext ? endNode.GetLeftExtensionNodesWithOrientation() : endNode.GetRightExtensionNodesWithOrientation();
                    // If this path diverges, we don't continue to follow it.
                    if (nextNodes.Count != 1)
                    {
                        finishedPaths.Add(path);
                        continue;
                    }
                    KeyValuePair<DeBruijnNode, bool> nextNode = nextNodes.First();
                    if (path.Nodes.Contains(nextNode.Key))
                    {
                        // Loop in path
                        //TODO: Not necessarily true, could overlap with itself but go out the other way
                        finishedPaths.Add(path);
                        continue;
                    }
                    else
                    {
                        // Update path orientation
                        path.IsSameOrientation = !(path.IsSameOrientation ^ nextNode.Value);
                        path.Nodes.Add(nextNode.Key);
                    }
                }
                maxDivergingPathLength++;
               
                /* Now to check for convergence, this is true if all paths can end with the same node
                   equivalent to all paths having the same node somewhere. 
                   TODO: Slow implementation is brute force N by all measure
                   first step would be to search only over the smallest possible path */
                var firstPathNodes = divergingPaths[0].Nodes;
                DeBruijnNode endingNode = null;
                for (int i = 1; i < firstPathNodes.Count;i++)
                {
                    var presentInAll = true;
                    var cur_node = firstPathNodes[i];
                    for (int k = 1; k <  divergingPaths.Count; k++) {
                        var c_path = divergingPaths[k];
                        if (!c_path.Nodes.Contains(cur_node))
                        {
                            presentInAll = false;
                            break;
                        }
                    }
                    if (presentInAll) 
                    {
                        endingNode = cur_node;
                        break;
                    }
                }
                // Paths have been extended. Check for convergence
                if (endingNode != null)
                {
                    DeBruijnPathList dpl = new DeBruijnPathList(divergingPaths.Count);
                    //If they have all converged, we now trim off any nodes at the end that didn't apply.
                    for (int i=0; i<divergingPaths.Count;i++)
                    {
                        var cur_path = divergingPaths[i];
                        DeBruijnPath dp;
                        if (endingNode != cur_path.Nodes.Last())
                        {
                            var indexOfEnd = cur_path.Nodes.IndexOf(endingNode);
                            dp = new DeBruijnPath(cur_path.Nodes.Take(indexOfEnd + 1));                            
                        }
                        else
                        {
                            dp =new DeBruijnPath(cur_path.Nodes);
                        }
                        dpl.AddPath(dp);
                    }

                    // Note: all paths have the same end node.
                    lock (redundantPaths)
                    {
                        // Redundant paths found
                        redundantPaths.Add(dpl);
                    }
                    return;
                }
            }
        }
    }
}
