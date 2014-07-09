using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;
using Bio.Algorithms.Alignment;
using Bio.Algorithms.Assembly.Graph;
using MitoDataAssembler;
using Bio;
using System.IO;

namespace MitoDataAssembler.Visualization
{

    public class GraphGenerator
    {
        DeBruijnGraph _graph;
        public List<MetaNode> MetaNodes = new List<MetaNode>();
        public DeBruijnGraph Graph { get { return _graph; } }
        public GraphGenerator(DeBruijnGraph assemblyGraph)
        {           
            this._graph = assemblyGraph;
            CreateMetaNodes();
            //verify all nodes visited
            if (_graph.GetNodes().Any(x => !x.IsVisited))
            {
                throw new Exception("Failed to visit all nodes!");
            }
        }
        public GraphGenerator(MitochondrialAssembly possibleAssembly) 
        {
            NodeCount = 0;
            this.MetaNodes = possibleAssembly.AllNodesInGraph;
        }

        /// <summary>
        /// Condense redundant paths down to simple paths
        /// </summary>
        /// <returns>List of simple paths.</returns>
        private void CreateMetaNodes()
        {
            _graph.SetNodeVisitState(false);
            //First step now, condense all nodes into "MetaNodes" that are linearly connected.
            //Note: Loop avoids stack overflow.
            foreach (DeBruijnNode node in _graph.GetNodes())
            {
                if (node.IsVisited) continue;
                else
                {
                    var megaNode = new MetaNode(node, _graph);
                    MetaNodes.Add(megaNode);                    
                }
            }
        }
  
        public static int NodeCount = 0;
        StreamWriter SW;
        private void OutNode(MetaNode node)
        {
            SW.Write(@"<node id='n" + node.NodeNumber.ToString() + @"'>
      <data key='V-AvgKmerCoverage'>" + node.AvgKmerCoverage.ToString() + @"</data>
      <data key='V-NodeCount'> " + node.ConstituentNodes.Count.ToString() + @"</data>
      <data key='V-Sequence'>" + node.Sequence + @"</data>
      <data key='V-LowRefPos'>" + node.Lowest_Reference_Position.ToString() + @"</data>
      <data key='V-HighRefPos'>" + node.Highest_Reference_Position.ToString() + @"</data>
    </node>
");
        }
        private void OutEdge(MetaNode n1, MetaNode n2,bool DifferentOrientation,int Right,uint weight)
        {
            
            SW.Write(@"<edge source='n" + n1.NodeNumber.ToString() + "' target='n" + n2.NodeNumber.ToString() + @"'>
            <data key='E-DifferentOrientation'>" + DifferentOrientation.ToString() + @"</data>
             <data key='E-Weight'>" + weight.ToString() + @"</data>
             <data key='E-GoingRight'>" + Right.ToString() + @"</data>
             </edge>
            ");
        }

        public void OutputDotGraph(string fname)
        {
            StreamWriter SW = new StreamWriter(fname);
            SW.WriteLine("digraph g {");
            foreach (MetaNode node in MetaNodes)
            {
                if (node.Lowest_Reference_Position != 0)
                {
                    SW.WriteLine("N" + node.NodeNumber.ToString() + " [label=\"" + node.Lowest_Reference_Position.ToString() + "-" + node.Highest_Reference_Position.ToString() + " -AvgCov="+node.AvgKmerCoverage.ToString()+ "\"];");
                }
                else
                {
                    SW.WriteLine("N" + node.NodeNumber.ToString() + " [label=\"" + node.Sequence + " -N="+node.AvgKmerCoverage.ToString()+"\"];");
                }
                }
            int edgeCount = 0;
            //make edges
            foreach (MetaNode node in MetaNodes)
            {

                foreach (var leftNode in node.GetLeftEdges().Where(x=>!x.IsInferiorEdge))
                {
                    edgeCount++;
                    string edge="N"+leftNode.FromNode.NodeNumber.ToString() + " -> "+"N"+leftNode.ToNode.NodeNumber.ToString()+";";
                    SW.WriteLine(edge);
                }
                //Remove all relatives
                foreach (var rightNode in node.GetRightEdges().Where(x=>!x.IsInferiorEdge))
                {

                    edgeCount++;
                    string edge = "N" + rightNode.FromNode.NodeNumber.ToString() + " -> " + "N" + rightNode.ToNode.NodeNumber.ToString() + ";";
                    SW.WriteLine(edge);
                }
                

            }
            SW.WriteLine("}");
            MitoPaintedAssembler.RaiseStatusEvent("\tWrote " + edgeCount.ToString() + " edges in dot file");
            SW.Close();
            
        }

        public void OutputGraph(string fname)
        {
            SW = new StreamWriter(fname);
            MitoPaintedAssembler.RaiseStatusEvent("\tBuilding Graph File");
            SW.Write(HEADER);
            //make verticies
            foreach (MetaNode node in MetaNodes)
            { OutNode(node);}
            int edgeCount = 0;
            //make edges
            foreach (MetaNode node in MetaNodes)
            {
                foreach (var leftNode in node.GetLeftEdges())
                {
                    edgeCount++;
                    OutEdge(node, leftNode.ToNode, leftNode.DifferentOrientation, 0, leftNode.Weight);
                }
                //Remove all relatives
                foreach (var rightNode in node.GetRightEdges())
                {
                    edgeCount++;
                    OutEdge(node, rightNode.ToNode, rightNode.DifferentOrientation, 100, rightNode.Weight);
                }
            }
            MitoPaintedAssembler.RaiseStatusEvent("\tWrote " + edgeCount.ToString() + " edges");
            SW.Write(FOOTER);
            SW.Close();
                    }

        public string HEADER = @"<?xml version='1.0' encoding='UTF-8'?>
<graphml xmlns='http://graphml.graphdrawing.org/xmlns'>
  <key id='V-AvgKmerCoverage' for='node' attr.name='AvgKmerCoverage' attr.type='double' />
  <key id='V-NodeCount' for='node' attr.name='NodeCount' attr.type='double' />
  <key id='V-Sequence' for='node' attr.name='Sequence' attr.type='string' />
  <key id='V-LowRefPos' for='node' attr.name='LowRefPos' attr.type='integer' />
  <key id='V-HighRefPos' for='node' attr.name='HighRefPos' attr.type='integer' />
  <key id='E-DifferentOrientation' for='edge' attr.name='DifferentOrientation' attr.type='boolean' />
  <key id='E-Weight' for='edge' attr.name='Weight' attr.type='double' />
  <key id='E-GoingRight' for='edge' attr.name='GoingRight' attr.type='double' />
  <graph edgedefault='directed'>
";

        public string FOOTER = @"  </graph>
</graphml>";

       
    }
}
