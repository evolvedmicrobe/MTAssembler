using System.Collections.Generic;

namespace HaploGrepSharp
{
    public class PhyloTreePath
    {
        private List<PhyloTreeNode> path = new List<PhyloTreeNode>();
        public PhyloTreePath(PhyloTreePath usedPath)
        {
            this.path.AddRange(usedPath.Nodes);
        }
        public PhyloTreePath()
        {
        }
        public virtual List<PhyloTreeNode> Nodes
        {
            get
            {
                return this.path;

            }
        }
        public virtual void add(PhyloTreeNode newNode)
        {
            this.path.Add(newNode);

        }
    }
}