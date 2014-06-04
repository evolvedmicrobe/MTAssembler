using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MitoDataAssembler.Visualization
{
    public class Arrow
    {
        public double sMag, sAng, eMag, eAng;
        internal MetaNode.Edge Edge;
        public Arrow(MetaNode.Edge edge)
        {
            this.Edge = edge;
            //need to set variables
            var toNode = edge.ToNode;
            var fromSlice = edge.FromNode.parentSlice;
            var toSlice = edge.ToNode.parentSlice;
            sAng = edge.GoingRight ? fromSlice.ActualRightGraphicalPosition.Value : fromSlice.ActualLeftGraphicalPosition.Value;
            eAng = edge.GoingRight ? toSlice.ActualLeftGraphicalPosition.Value : toSlice.ActualRightGraphicalPosition.Value;
            sMag = fromSlice.ActualHeightMidpoint;
            eMag = toSlice.ActualHeightMidpoint;
        }
        
        
    }
}
