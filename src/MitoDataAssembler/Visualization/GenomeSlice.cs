﻿using MitoDataAssembler;
using MitoDataAssembler.Visualization;

namespace MitoDataAssembler.Visualization
{
    /// <summary>
    /// Represent a slice of a section of the genome to be plotted, basically a metanode with enhance capabilities
    /// </summary>
    public class GenomeSlice
    {
        public string Label;

        private double? _suggestedRightSidePosition;
        private double? _suggestedLeftSidePosition;

        /// <summary>
        /// Returns the suggested position for the right side, actual position will be the center of the suggested positions
        /// plus the sequence length
        /// </summary>
        public double? SuggestedRightSidePosition
        {
            get
            {
                double? right = ActualRightGraphicalPosition;
                if (right.HasValue)
                {
                    return right.Value;
                }
                else
                {
                    return _suggestedRightSidePosition;
                }
            }
            set
            {
                _suggestedRightSidePosition = value;
            }

        }
        public double? SuggestedLeftSidePosition
        {
            get
            {
                double? left = ActualLeftGraphicalPosition;
                if (left.HasValue)
                {
                    return left.Value;
                }
                else
                {
                    return _suggestedLeftSidePosition;
                }
            }
            set
            {
                _suggestedLeftSidePosition = value;
            }
        }
        public double? ActualLeftGraphicalPosition
        {
            get
            {
                if (_suggestedLeftSidePosition.HasValue && _suggestedRightSidePosition.HasValue)
                {
                    double center = calculateCenter();
                    double midlength = pnode.LengthOfNode / 2.0;
                    var left = center - midlength;
                    if (left > Utils.CircularGenomeCaseHandlers.MT_Genome_Length)
                        return left - Utils.CircularGenomeCaseHandlers.MT_Genome_Length;
                    else if (left < 0)
                        return Utils.CircularGenomeCaseHandlers.MT_Genome_Length + left;
                    else
                        return left;
                }
                else
                    return null;
            }
        }
        public double? ActualRightGraphicalPosition
        {
            get
            {
                if (_suggestedLeftSidePosition.HasValue && _suggestedRightSidePosition.HasValue)
                {
                    //are we flipped? 

                    double center = calculateCenter();
                    double midlength = this.pnode.LengthOfNode / 2.0;
                    var right = center + midlength;
                    if (right > Utils.CircularGenomeCaseHandlers.MT_Genome_Length)
                        return right - Utils.CircularGenomeCaseHandlers.MT_Genome_Length;
                    else
                        return right;
                }
                else
                    return null;
            }
        }

        /// <summary>
        /// Does the segment span the end and loop back? Such that left > right
        /// </summary>
        /// <returns></returns>
        public bool SegmentLoopsAround()
        {
            return Utils.CircularGenomeCaseHandlers.SegmentLoopsAround(_suggestedLeftSidePosition.Value, _suggestedRightSidePosition.Value,Node.LengthOfNode );
        }
        public bool GraphSegmentsLoopAround()
        {
            return Utils.CircularGenomeCaseHandlers.SegmentLoopsAround(this.ActualLeftGraphicalPosition.Value, this.ActualRightGraphicalPosition.Value, Node.LengthOfNode);
        }

        private double? pActualCenter;
        public double ActualHeightMidpoint
        {
            get {
                return MagnitudeLow + Height * .5; 
            }
        }
        public double CalculateRealizedCenter()
        {
            if (!pActualCenter.HasValue)
            {
                if (!SegmentLoopsAround())
                {
                    pActualCenter = (ActualLeftGraphicalPosition + ActualRightGraphicalPosition) / 2.0;
                }
                else
                {
                    var temp = (Utils.CircularGenomeCaseHandlers.MT_Genome_Length - ActualLeftGraphicalPosition) + ActualRightGraphicalPosition;
                    temp = temp / 2.0;
                    temp = temp + ActualLeftGraphicalPosition;
                    if (temp > Utils.CircularGenomeCaseHandlers.MT_Genome_Length)
                        temp = temp - Utils.CircularGenomeCaseHandlers.MT_Genome_Length;
                    pActualCenter = temp;
                }
            }
            return pActualCenter.Value;
        }
        private double calculateCenter()
        {
            if (!SegmentLoopsAround())
            {
                return (_suggestedLeftSidePosition.Value + _suggestedRightSidePosition.Value) / 2.0;
            }
            else
            {
                //Otherwise, might be in the middle somewhere
                bool leftLow = _suggestedLeftSidePosition.Value > Utils.CircularGenomeCaseHandlers.MT_Genome_Length / 2.0;
                double left = leftLow ? _suggestedLeftSidePosition.Value : _suggestedRightSidePosition.Value;
                double right = leftLow ? _suggestedRightSidePosition.Value : _suggestedLeftSidePosition.Value;
                right = right + Utils.CircularGenomeCaseHandlers.MT_Genome_Length;
                //should always be the length of the guy
                return (left + right) / 2.0;
            }
        }
     
            public double Height;
            public bool InMainAssembly;
            public double MagnitudeLow;
            public double MagnitudeHigh { get { return MagnitudeLow + Height; } }
     
        /// <summary>
        /// Gets or sets Fill.
        /// </summary>
        public string Fill
        {
            get
            {
                if (SuggestedRightSidePosition.Value < SuggestedLeftSidePosition.Value &&
                    !SegmentLoopsAround())
                {
                    return "blue";
                }
                else
                {
                    return "red";
                }
            }
        }
        private MetaNode pnode;
        public MetaNode Node {get{return pnode;}}
        /// <summary>
        /// Initializes a new instance of the <see cref="PieSlice"/> class.
        /// </summary>
        /// <param name="label">
        /// The label.
        /// </param>
        /// <param name="value">
        /// The value.
        /// </param>
        /// <param name="fill">
        /// The fill.
        /// </param>
        public GenomeSlice(MetaNode node)
        {
            this.Label = "N-"+node.NodeNumber.ToString();
            this.pnode=node;
            node.parentSlice = this; 
        }
      
        
    }
}