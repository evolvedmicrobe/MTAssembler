using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;

namespace HaploGrepSharp
{
    public class SampleFile
    {
        internal Dictionary<string, TestSample> testSamples = new Dictionary<string, TestSample>();
        public SampleFile(List<string> sampleLines)
        {
            int lineIndex = 1;
            foreach (string currentLine in sampleLines)
            {
                TestSample newSample;
                try
                {
                    newSample = new TestSample(currentLine);
                }
                catch (Exception e)
                {
                    throw new HaploGrepException("Could not parse line " + lineIndex.ToString(), e);
                }

                if (this.testSamples.ContainsKey(newSample.SampleID))
                { throw new HaploGrepException("Two samples appear to have the same name: " + newSample.SampleID + " second one found at line: " + lineIndex.ToString()); }
                this.testSamples[newSample.SampleID] = newSample;
                lineIndex++;
            }
        }

        private SampleRange tryDetectRange(string line)
        {
            throw new Exception("Not implemented");
            //if (line.Contains("#!"))
            //{
            //    SampleRange range = new SampleRange();
            //    StringTokenizer rangeTokenizer = new StringTokenizer(line, " ");
            //    rangeTokenizer.nextToken();

            //    while (rangeTokenizer.hasMoreElements())
            //    {
            //        string rangeToken = rangeTokenizer.nextToken().Trim();
            //        if (rangeToken.Contains("-"))
            //        {
            //            range.addCustomRange((int)Convert.ToInt32(rangeToken.Substring(0, rangeToken.IndexOf("-"))), (int)Convert.ToInt32(rangeToken.Substring(rangeToken.IndexOf("-") + 1, rangeToken.Length - (rangeToken.IndexOf("-") + 1))));


            //        }
            //        else
            //        {
            //            range.addCustomRange((int)Convert.ToInt32(rangeToken), (int)Convert.ToInt32(rangeToken));
            //        }

            //    } return range;

            //}
            //return null;

        }

        public virtual List<TestSample> TestSamples
        {
            get
            {
                return new List<TestSample>(this.testSamples.Values);

            }
        }
        public virtual string ToString()
        {
            string result = "";

            foreach (TestSample currenTestSample in this.testSamples.Values)
            {
                result = result + currenTestSample.ToString() + "\n";

            }
            return result;

        }
    }
}