using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
namespace HaploGrepSharp
{
    public class TestSample : IComparable<TestSample>
    {
        private Haplogroup predefiniedHaplogroup;
        private Haplogroup recognizedHaplogroup;
        private Sample sample;
        private string testSampleID = "Unknown";
        private SampleRange sampleRange = null;
        private string state = "n/a";
        private double resultQuality = 0.0D;
        public TestSample(string sampleID, Haplogroup predefiniedHaplogroup, Sample sample, SampleRange sampleRange, string state)
        {
            this.testSampleID = sampleID;
            this.sample = sample;
            this.predefiniedHaplogroup = predefiniedHaplogroup;
            this.sampleRange = sampleRange;
            this.state = state;
        }
        public TestSample(string currentLine)
        {
            string[] tokens = currentLine.Split(new char[] { '\t' }, StringSplitOptions.RemoveEmptyEntries);
            if (tokens.Length < 4) {
                throw new HaploGrepException("Not enough columns in .hsd file");
            }
            this.testSampleID = tokens[0].Trim();
            if (tokens[1].Length == 0) {
                throw new ArgumentException("No range specified, range given was: " + tokens[1]);
            }
            tokens[1] = tokens[1].Replace("\"", "");
            this.sampleRange = new SampleRange(tokens[1]);
            if ((tokens[2].Equals("?")) || (tokens[2].Equals("SEQ"))) {
                this.predefiniedHaplogroup = new Haplogroup("");
            }
            else {
                this.predefiniedHaplogroup = new Haplogroup(tokens[2]);
            }
            StringBuilder sampleString = new StringBuilder();
            for (int i = 3; i < tokens.Length; i++)
            {
                sampleString.Append(tokens[i] + " ");
            }
            try
            {
                this.sample = new Sample(sampleString.ToString(), 0);
            }
            catch (Exception e)
            {
                throw e;
            }
        }
        public Haplogroup PredefiniedHaplogroup
        {
            get
            {
                return this.predefiniedHaplogroup;
            }
            set
            {
                this.predefiniedHaplogroup = value;
            }
        }
        public List<Polymorphism> Polymorphismn
        {
            get
            {
                return this.sample.sample;
            }
        }
        public SampleRange SampleRanges
        {
            get
            {
                return this.sampleRange;
            }
        }
        public Sample Sample
        {
            get
            {
                return this.sample;
            }
        }
        public override string ToString()
        {
            string result = this.testSampleID + "\t" + this.sampleRange + "\t" + this.predefiniedHaplogroup + "\t";
            foreach (Polymorphism currentPoly in this.sample.sample)
            {
                result = result + currentPoly.ToString() + " ";
            }
            return result;
        }
        public string SampleID
        {
            get {
                return this.testSampleID;
            }
        }
        public string Status
        {
            get{
                return this.state;
            }
        }
        public string State
        {
            set
            {
                this.state = value;
            }
        }
        public Haplogroup RecognizedHaplogroup
        {
            get
            {
                return this.recognizedHaplogroup;
            }
            set
            {
                this.recognizedHaplogroup = value;
            }
        }

        public double ResultQuality
        {
            set
            {
                this.resultQuality = value;
            }
            get
            {
                return this.resultQuality;
            }
        }
        public List<Polymorphism> PolyNotinRange
        {
            get
            {
                List<Polymorphism> notInRangePolys = new List<Polymorphism>();
                foreach (Polymorphism currentPoly in this.sample.Polymorphismn)
                {
                    if (!this.sampleRange.contains(currentPoly))
                    {
                        notInRangePolys.Add(currentPoly);
                    }
                } return notInRangePolys;
            }
        }
        public int CompareTo(TestSample o)
        {
            if (SampleID.CompareTo(o.SampleID) < 0)
            {
                return -1;
            }
            if (SampleID.CompareTo(o.SampleID) > 0)
            {
                return 1;
            } return 0;
        }
    }
}