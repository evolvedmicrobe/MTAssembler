using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using Bio;
using Bio.IO.BAM;
using Bio.IO.SAM;
using System.IO;
using Bio.Algorithms.Assembly.Padena;
using RDotNet;

namespace MitoDataAssembler.PairedEnd
{
	/// <summary>
	/// This class extracts paired end mtDNA reads from a BAM file and looks for peaks in the template length that are indicative of 
	/// a deletion by looking for peaks.
	/// </summary>
	public class PairedEndPeakFinder
	{
		/// <summary>
		/// If less than this number of read pairs are observed in the search range, skip the match.
		/// </summary>
		public const double MINIMUM_REQUIRED__READS_FOR_SEARCH=100;

		/// <summary>
		/// Local peaks are called at this level if there occurence is less than this value
		/// across all tested sites.
		/// </summary>
		private double FamilyAlphaValue = 0.05;

		/// <summary>
		/// The poisson value used to test for the signficance level required for a peak.
		/// </summary>
		/// <value>The poisson lambda.</value>
		[OutputAttribute]
		public double PoissonLambda { get; private set; }

		/// <summary>
		/// Smoothing window size used for local minima detection.
		/// </summary>
		[OutputAttribute]
		public int PeakEdgeWindowSize {get;private set;}

        /// <summary>
        /// The value of the window size used to scan for deletions, the average number of items
        /// in this window will be checked.  Should be odd
        /// </summary>
        [OutputAttribute]
        public int ScanningWindowSize {get;private set;} 

		/// <summary>
		/// The file name used
		/// </summary>
		[OutputAttribute]
		public string FileName{ get; private set; }

		/// <summary>
		/// The mtDNA used for querying the BAM
		/// </summary>
		/// <value>The name of the mt DNA.</value>
		[OutputAttribute]
		public string mtDNAName{ get; private set; }

		/// <summary>
		/// How many total reads?
		/// </summary>
		[OutputAttribute]
		public int TotalMTReads{ get; private set; }

		/// <summary>
		/// Count that had insert size over 0
		/// </summary>
		[OutputAttribute]
		public int CountZero{ get; private set; }

		/// <summary>
		/// Those over the CRS in length
		/// </summary>
		[OutputAttribute]
		public int CountOverMax{ get; private set; }

		/// <summary>
		/// The lowest value from the tail to look for deletions in.  This should be the tail of the 
		/// distribution that corresponds to the "real" deletion sizes.
		/// </summary>
		[OutputAttribute]
		public int? LowEndOfSearchRange { get; private set; }

		/// <summary>
		/// Gets the high end of search range.
		/// </summary>
		/// <value>The high end of search range.</value>
		[OutputAttribute]
		public int? HighEndOfSearchRange{ get; private set; }

        [OutputAttribute]
        public string ReasonForFailure { get; private set; }

        [OutputAttribute]
        public PossibleDeletionCollection DeletionSizesFromPeakFinding {get; private set;}

		/// <summary>
		/// A value determines if the search could feasably be done.
		/// </summary>
		/// <value><c>true</c> if search was possible; otherwise, <c>false</c>.</value>
		[OutputAttribute]
		public bool SearchWasPossible { get; private set; }

		/// <summary>
		/// Size of the most frequently observed template length
		/// </summary>
		/// <value>The max inset size observed.</value>
		[OutputAttribute]
		public int MostFrequentInsetSizeObserved{ get; private set; }

		/// <summary>
		/// Number of times the maximally observed insert size appeared.
		/// </summary>
		/// <value>The max value.</value>
		[OutputAttribute]
		public int MaximumObservedCounts { get; private set; }

        /// <summary>
        /// Number of possible deletions found
        /// </summary>
        [OutputAttribute]
        public int PossibleDeletionCount { get { return DeletionSizesFromPeakFinding.Count; } }

        /// <summary>
        /// Cutoff used to find peaks, anything above or equal to this value counts.
        /// </summary>
        [OutputAttribute]
        public double PoissonCoverageCutoffCount { get; private set; }

        /// <summary>
        /// The prefix to use in front of output files
        /// </summary>
        [OutputAttribute]
        public string OutputFilePrefix { get; private set; }

		BAMParser bp;
        RInterface rInt;
		/// <summary>
		/// Position i is the count of reads with insert size i
		/// </summary>
		double[] countsOfSize;

		public PairedEndPeakFinder (string BAMFile,string prefix,string mtDNAName="MT")
		{
         //   ReasonForFailure = "NOT SET";
            PeakEdgeWindowSize = 3;
            ScanningWindowSize = 11;
			if(String.IsNullOrEmpty(BAMFile) || String.IsNullOrEmpty(mtDNAName)  || String.IsNullOrEmpty(prefix))
			{
				throw new ArgumentException("Tried to load the paired peak finder without the a BAM file, output file Prefix or mt DNA chromosome name");
			}
			if(!File.Exists(BAMFile))
			{
				throw new IOException("Could not find BAM file: "+BAMFile);
			}
            rInt = new RInterface();
			FileName = BAMFile;
			this.mtDNAName = mtDNAName;
            this.OutputFilePrefix = prefix;
			bp = new BAMParser ();
			//create count array
			countsOfSize = new double[StaticResources.CRS_LENGTH];

		}
        
		/// <summary>
		/// Read the file and scan for deletions.
		/// </summary>
		public void FindDeletionPeaks()
		{
            DeletionSizesFromPeakFinding = new PossibleDeletionCollection();
            SearchWasPossible = true;
			//Load the distribution
			loadDistribution ();
			//set the range to search (low to high along chromosome)
			setSearchRange ();
            if (SearchWasPossible)
            {
                double sumInRange = GetRangeToEvaluate().Sum();
                if (sumInRange < MINIMUM_REQUIRED__READS_FOR_SEARCH)
                {
                    SearchWasPossible = false;
                    ReasonForFailure = "Number of reads in search range was less than minimum allowed.  Found: " + sumInRange.ToString() + " Needed: " + MINIMUM_REQUIRED__READS_FOR_SEARCH.ToString();
                }
                else
                {
                    //set poisson threshold out 
                    setPoissonThreshold();
                    Console.WriteLine("Attempting deletion finding");
                    DeletionSizesFromPeakFinding.AddRange(ScanForDeletions());
                    Console.WriteLine("Ending deletion finding");
                }
                //make a plot with the values
                makeTotalCoveragePlotOverSearchRange();
				makeTotalCoveragePlot ();
                //also do it for window averaged
                makeSlidingWindowOccurencePlot();
                    
            }
		}
   

        /// <summary>
        /// Scan sliding windows, merge all overlapping windows that contain a peak, then return these intervals.
        /// </summary>
        /// <returns></returns>
        IEnumerable<PossibleDeletion> ScanForDeletions()
        {
            /* Calculate amount to merge non-overlapping windows, 
             * use the width of the observed template size distribution
             * is used as the amount
             */
            var MergeDistance = LowEndOfSearchRange.Value - MostFrequentInsetSizeObserved;

            //get list of all windows with average coverage above the threshold
            var windows = GetAverageInSlidingCoverageWindows().Where(z => z.Value > PoissonCoverageCutoffCount).ToList();
            PossibleDeletion toRet = null;
            int lastPos = -10;
            foreach (var win in windows)
            {
                //we have a peak, let's find the maximum value (assuming window size in max is 
                if (toRet == null)
                {
                    toRet = new PossibleDeletion() { Start = win.Key };
                    lastPos = win.Key;
                }
                else if ((win.Key - lastPos) < MergeDistance)
                {
                    toRet.End = win.Key;
                    lastPos = win.Key;
                }
                else
                {
                    if (toRet.End == 0) { toRet.End = toRet.Start; }
                    yield return toRet;
                    lastPos = win.Key;
                    toRet = new PossibleDeletion() { Start = win.Key };
                }

            }
            if (toRet != null)
            {
                if (toRet.End == 0) { toRet.End = toRet.Start; }
                yield return toRet;
            }
        }
        public IEnumerable<KeyValuePair<int, double>> GetAverageInSlidingCoverageWindows()
        {
            int curPos = LowEndOfSearchRange.Value;
            int high = HighEndOfSearchRange.Value;
            var avg = countsOfSize.Skip(curPos).Take(ScanningWindowSize).Sum();
            var delta = ScanningWindowSize/2;
            var dwin = (double)ScanningWindowSize;
            yield return new KeyValuePair<int,double>(curPos+delta,avg/dwin);
            curPos++;
            while(curPos<(HighEndOfSearchRange-ScanningWindowSize))
            {
                //remove the last value, add the next one
                avg = avg - countsOfSize[curPos - 1] + countsOfSize[curPos + ScanningWindowSize-1];
                yield return new KeyValuePair<int, double>(curPos + delta, avg / dwin);
                curPos++;
            }
        }
        /// <summary>
        /// Loads the distribution of insert sizes, creates an array of [counts]
        /// </summary>
        void loadDistribution()
		{
			//load all the data
			var producer=bp.ParseRangeAsEnumerableSequences (FileName, mtDNAName, 0, Int32.MaxValue);
			//Go through and make an array of counts
			foreach (var s in producer) {
				TotalMTReads++;
				if ((s.Flag & SAMFlags.PairedRead) != SAMFlags.PairedRead) {
					continue;
				} else if (s.ISize == 0) {
					CountZero++;
				} else {
					int i = s.ISize;
					//those less than 0 should have their mate contain a value, can continue
					if (i < 0) {
						continue;
					} else if (i > (StaticResources.CRS_LENGTH - 1)) {
						CountOverMax++;
					} else {
						countsOfSize [i]+=1.0;
					}
				}
			}

		}
		/// <summary>
		/// Scan the distribution to define the range over which we will be searching.
		/// </summary>
		private void setSearchRange()
		{
			double max = countsOfSize.Max ();
			int maxIndex;
			//find the max value, and the first local minimum after it, assumes function is monotonic
			for (maxIndex=0; maxIndex < countsOfSize.Length; maxIndex++) {
				if (countsOfSize[maxIndex] == max) {
					break;
				}
			}
			//proceed forward to find first local minima, go in windows of three to smooth somewhat.  This is the end of the
			//bulk of the insert size distribution.
			double lastVal = max;
			for (int i=maxIndex+PeakEdgeWindowSize; i < (countsOfSize.Length-PeakEdgeWindowSize); i+=PeakEdgeWindowSize) {
				double average = 0;
				for (int j=0; j<PeakEdgeWindowSize; j++) {
					average += countsOfSize [i+j];
				}
				average = average / PeakEdgeWindowSize;
				if (average >= lastVal || average == 0) {
					LowEndOfSearchRange = i + 1;
					break;
				}				
			}
			if (!LowEndOfSearchRange.HasValue) {
				SearchWasPossible = false;
                ReasonForFailure = "Low end of search range was never set";
			} else {
				//set high end of the search by simply reflecting this value over the end;
				HighEndOfSearchRange = countsOfSize.Length - LowEndOfSearchRange;
				MostFrequentInsetSizeObserved = maxIndex;
				MaximumObservedCounts = (int)max;
			}
		}

        /// <summary>
        /// Get just the interval of the genome to be covered by the mtDNA.
        /// </summary>
        /// <returns></returns>
        private List<double> GetRangeToEvaluate()
        {
            return countsOfSize.Skip(LowEndOfSearchRange.Value).Take(HighEndOfSearchRange.Value - LowEndOfSearchRange.Value + 1).ToList();
        }

        /// <summary>
        /// Plot the intervals and densities as we move along.
        /// </summary>
        private void makeTotalCoveragePlotOverSearchRange()
        {
            var counts = GetRangeToEvaluate();
            var engine = rInt.CurrentEngine;
            lock (engine)
            {

                var yvals = engine.CreateNumericVector(counts);
                var xvals = engine.CreateNumericVector(Enumerable.Range(LowEndOfSearchRange.Value, counts.Count).Select(x => (double)x));
                string fname = OutputFilePrefix + "SearchRange.pdf";
                var scutoff = PoissonCoverageCutoffCount.ToString();
                string line = "lines(c(0,16569),c(" + scutoff + "," + scutoff + "),col=\"red\",lwd=3)";
                rInt.PlotPDF(xvals, yvals, fname, "Template Length Distribution", "Length", "Counts", new List<string>() { line });
            }
        }

		private void makeTotalCoveragePlot()
		{
			var lowVal = LowEndOfSearchRange.Value;
			var highVal = HighEndOfSearchRange.Value;
			var y_high = (countsOfSize.Max () * 5).ToString(); 
			var engine = rInt.CurrentEngine;
			lock (engine)
			{
				var xvals = engine.CreateNumericVector(Enumerable.Range(1,countsOfSize.Length).Select(x=>(double)x));
				var yvals = engine.CreateNumericVector(countsOfSize);
				string fname = OutputFilePrefix + "FullRange.pdf";
				var scutoff = PoissonCoverageCutoffCount.ToString();
				List<string> addCmds = new List<string> (3);
				string line = "lines(c(0,16569),c(" + scutoff + "," + scutoff + "),col=\"red\",lwd=3)";
				addCmds.Add (line);
				line = "lines(rep("+lowVal.ToString()+",2),c(0, " + y_high + "),col=\"red\",lwd=3)";
				addCmds.Add (line);
				line = "lines(rep("+highVal.ToString()+",2),c(0, " + y_high + "),col=\"red\",lwd=3)";
				addCmds.Add (line);
				rInt.PlotPDF(xvals,yvals, fname,"Template Length Distribution", "Template Length", "Counts of Reads", addCmds);

			}
		}
        private void makeSlidingWindowOccurencePlot()
        {
            MitoPaintedAssembler.RaiseStatusEvent("\tStarting sliding window plot");
            var list = GetAverageInSlidingCoverageWindows().ToList();
            double[] xvals = list.Select(x => (double)x.Key).ToArray();
            double[] yvals = list.Select(x => x.Value).ToArray();
            var scutoff = PoissonCoverageCutoffCount.ToString();
            string line = "lines(c(0,16569),c(" + scutoff + "," + scutoff + "),col=\"red\",lwd=3)";
            rInt.PlotPDF(xvals, yvals, OutputFilePrefix + "Windowed" + this.ScanningWindowSize.ToString() + ".pdf", "Windowed occurence", "Insert Size", "Window Average", new List<string>() { line });
            MitoPaintedAssembler.RaiseStatusEvent("\tFinished sliding window plot");

        }
        /// <summary>
		/// Sets the threshold by assuming a poisson distribution over the "normal" values whose distribution is
		/// assumed to be a poisson that is estimated robustly using the median of counts as values.
		/// </summary>
		private void setPoissonThreshold()
		{
            //grab a list of counts in this range
            var countCopy = GetRangeToEvaluate();
			countCopy.Sort ();
			//grab median value
			var median = countCopy [countCopy.Count / 2];
            //if median 0, set value to average of zero and 1 values in the list
            if (median < .0001)
            {
                var zeroOnes = countCopy.Where(x => x < 1.25).ToList();
                median = zeroOnes.Sum() / zeroOnes.Count;
            }
            var length = countCopy.Count;
            //find p-value cutoff for highest value likely seen amongst length multiple draws of a poisson with that median value
			var cdfCutoff = Math.Log (1-FamilyAlphaValue)/(double)length;
            //qpois(cdf,1,log.p=TRUE)

            var engine = rInt.CurrentEngine;
            lock (engine)
            {
                PoissonLambda = median;
                var lambda = engine.CreateNumericVector(new double[] { median, cdfCutoff });
                engine.SetSymbol("lambda", lambda);
                var cutoff = engine.Evaluate("qpois(lambda[2],lambda[1],log.p=TRUE)").AsNumeric().First();
                System.Diagnostics.Debug.WriteLine("Lambda: " + median.ToString());
                System.Diagnostics.Debug.WriteLine("Cutoff: " + cutoff.ToString());
                System.Diagnostics.Debug.WriteLine("QPois: " + cdfCutoff.ToString());
                PoissonCoverageCutoffCount = cutoff;
            }
		}  
        
        




	}
}

