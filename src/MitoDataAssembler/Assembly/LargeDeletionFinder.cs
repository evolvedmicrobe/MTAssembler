using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Bio.Algorithms.Assembly.Graph;
using Bio;
using System.Collections.Concurrent;
using MitoDataAssembler.Visualization;
using MitoDataAssembler.Extensions;
using System.Diagnostics;
using System.IO;

namespace MitoDataAssembler
{


    /// <summary>
    /// Given a circular genome and a painted deBruijin Graph, attempts to find places in the genome that are indicative of a large deletion
    /// that is where a path in the graph joins to sections that are too far apart or too short. 
    /// 
    /// Somewhat inefficient algorithm at present which should be fine for mitochondria data, it works by following every painted path from 
    /// a k-mer and reporting the Path if the difference is  > the threshold.  If so, it reports this junction +/- the window around it.
    /// 
    /// </summary>
    class LargeDeletionFinder
    {
        #region Fields, Constructor, Properties

        /// <summary>
        /// Holds reference to assembler graph.
        /// </summary>
        public static DeBruijnGraph graph;

        public List<PossibleAssembly> PossibleDeletionPaths = new List<PossibleAssembly>();
        public List<DeletionReport> DeletionReports;
        public static MitochondrialAssembly putativeAssembly;
        //public int MinUnExplainedDistance = 40;

        /// <summary>
        /// Initializes a new instance of the deletion finder, has the 
        /// </summary>
        /// <param name="lengthToReport"></param>
        public LargeDeletionFinder()
        {
        }

        /// <summary>
        /// Gets the name of the algorithm.
        /// </summary>
        public string Name
        {
            get { return "Mito Deletion Finder"; }
        }

        ///// <summary>
        ///// Gets or sets the minimum distance between two painted nodes before it is reported
        ///// </summary>
        //public int MinimumDistanceThreshold
        //{
        //    get;
        //    set;
        //}
        #endregion
        private int KmerLength;
        public List<DeletionReport> FindAllDeletions(MitochondrialAssembly assembly)
        {
            //look for sections weh

            return null;
        }
        public List<DeletionReport> FindAllDeletions(DeBruijnGraph graph, MitochondrialAssembly assembly)
        {
            LargeDeletionFinder.graph = graph;
            KmerLength = graph.KmerLength;
            //set all edges in the graph to not be visited
            graph.GetNodes().AsParallel().ForAll(x => x.ResetVisitState());
            foreach (DeBruijnNode node in graph.GetNodes())
            {
                //starting from any unused edges in the network, make any/all paths one can
				try
				{
                    PossibleDeletionPaths.AddRange(ExtendFromStartNode(node));
				}
				catch(Exception thrown) {
					Console.WriteLine (thrown.Message);
				}
            }
            DeletionReport.CompleteAssembly = assembly;
            DeletionReports = PossibleDeletionPaths.Select(x => new DeletionReport(x)).ToList();
            return DeletionReports;
        }

        public IEnumerable<PossibleAssembly> ExtendFromStartNode(DeBruijnNode start)
        {
            //TODO: I believe this handles figure 8s and palindromes just fine, should verify though.

            //First go Right
            var rightNeighbors = start.GetRightExtensionNodesWithOrientationMarkingEdgeAsVisited(false);
            List<PossibleAssembly> rights = new List<PossibleAssembly>();
            foreach (var direction in rightNeighbors)
            {
                PossibleAssembly pa = new PossibleAssembly(start, true);
                rights.AddRange(ExtendChain(pa, direction.Key, true, direction.Value));
            }
            List<PossibleAssembly> lefts = new List<PossibleAssembly>();
            var leftNeighbors = start.GetLeftExtensionNodesWithOrientationMarkingEdgeAsVisited(false);
            foreach (var direction in leftNeighbors)
            {
                PossibleAssembly pa = new PossibleAssembly(start, false);
                lefts.AddRange(ExtendChain(pa, direction.Key, false, direction.Value));
            }
            //Now to combine a left and right chain
            if (lefts.Count > 0 && rights.Count > 0)
            {
                foreach (var right in rights)
                {
                    foreach (var left in lefts)
                    {
                        yield return new PossibleAssembly(left, right);
                    }
                }
            }
            else if (lefts.Count > 0)
            {
                foreach (var left in lefts)
                {
                    yield return left;
                }
            }
            else if (rights.Count > 0)
            {
                foreach (var right in rights)
                {
                    yield return right;
                }
            }
        }
        private IEnumerable<PossibleAssembly> ExtendChain(PossibleAssembly currentPath, DeBruijnNode nextNeighbor, bool goingRight, bool sameOrientation)
        {   
            byte nextSymbol = MetaNode.GetNextSymbol(nextNeighbor, KmerLength, !goingRight);
            currentPath.Add(nextNeighbor,nextSymbol);
            nextNeighbor.IsVisited = true;
            bool nextRight = !goingRight ^ sameOrientation;
            List<KeyValuePair<DeBruijnNode, bool>> nextNodes = nextRight ? nextNeighbor.GetRightExtensionNodesWithOrientationMarkingEdgeAsVisited() :
            nextNeighbor.GetLeftExtensionNodesWithOrientationMarkingEdgeAsVisited();
            DeBruijnNode next;
            //DeBruijnNode last = currentPath.constituentNodes[currentPath.constituentNodes.Count-1];
            //DeBruijnNode first=currentPath.constituentNodes[0];
            while (nextNodes.Count == 1)
            {
                var nextSet = nextNodes.First();
                next = nextSet.Key;
                sameOrientation = nextSet.Value;
                nextRight = (!nextRight) ^ sameOrientation;
                nextSymbol = MetaNode.GetNextSymbol(next, KmerLength, !nextRight);
                //now check if we are in a circle or a loop at the end, these are very annoying situtations, basic criteria, can't leave
                //the same node the same way twice
                if (next.IsVisited && currentPath.constituentNodes.Contains(next))
                {

                    //okay, if we are equal to the first node or the last node, we can't leave or return the same way we came, otherwise we are done.
                    var excludedNextNodes = currentPath.GetPreviousWaysNodeWasLeft(next);
                    //how many neighbors dow we have in this group?
                    var temp = nextRight ? next.GetRightExtensionNodesWithOrientationMarkingEdgeAsVisited() : next.GetLeftExtensionNodesWithOrientationMarkingEdgeAsVisited();
                    temp = temp.Where(x => !excludedNextNodes.Contains(x.Key)).ToList();
                    //only one way to go
                    if (temp.Count == 1)
                    {
                        nextNodes = temp;
                        //currentPath.contigSequence.Add(nextSymbol);
                        currentPath.Add(next,nextSymbol);
                        next.IsVisited = true;//flag not actually used though 
                    }
                    else if (temp.Count == 0)//done
                    {
                        if (currentPath.constituentNodes[0] == next)
                        { currentPath.CircularLoop = true; }
                        yield return currentPath;
                        //nextNodes.Clear();//we are done
                        yield break;
                    }
                    else //Extend path using all feasible options, then continue.
                    {
                        foreach (var neighbor in temp)
                        {
                            foreach (var v in ExtendChain(currentPath.Clone(), neighbor.Key, nextRight, neighbor.Value))
                            {
                                yield return v;
                            }
                        }
                        //nextNodes.Clear();//done
                        yield break;
                    }
                }
                else
                {
                    //currentPath.contigSequence.Add(nextSymbol);
                    currentPath.Add(next,nextSymbol);
                    next.IsVisited = true;//flag not actually used though 
                    nextNodes = nextRight ? next.GetRightExtensionNodesWithOrientationMarkingEdgeAsVisited() : next.GetLeftExtensionNodesWithOrientationMarkingEdgeAsVisited();
                }
            }
            //If we have more than one node remaining, have to kick it off.
            if (nextNodes.Count > 1)
            {
                foreach (var neighbor in nextNodes)
                {
                    foreach (var v in ExtendChain(currentPath.Clone(), neighbor.Key, nextRight, neighbor.Value))
                    {
                        yield return v;
                    }
                }
            }
            if (nextNodes.Count == 0)
            {
                yield return currentPath;
            }
        }
        public void OutputReport(string name)
        {
            StreamWriter sw = new StreamWriter(name);
            sw.WriteLine(DeletionReport.DeletionReportHeaderLine());
            foreach (string str in
            DeletionReports.SelectMany(x => x.DeletionReportDataLines()))
            { sw.WriteLine(str); }
            sw.Close();
        }
    }
}

