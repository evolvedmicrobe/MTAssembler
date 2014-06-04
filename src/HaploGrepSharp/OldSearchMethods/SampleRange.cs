using System;
using System.Collections.Generic;
namespace HaploGrepSharp
{
    //TODO: Values should not be hard coded
    public class SampleRange
    {
        List<int> starts = new List<int>();
        List<int> ends = new List<int>();
        public List<int> Starts
        {
            get
            {
                return this.starts;
            }
        }
        public List<int> Ends
        {
            get
            {
                return this.ends;
            }
        }
        public SampleRange()
        {
        }
        public SampleRange(SampleRange rangeToCopy)
        {
            this.starts.AddRange(rangeToCopy.starts);
            this.ends.AddRange(rangeToCopy.ends);
        }
        public SampleRange(string rangesToParse)
        {
            //Strings are e.g. "16024-576;2092;3552;4071;4491;4833;4883;8414;8473;9090;9824;10397;10400;11959;11969;12372;12771;13563;14502;14569;15487"
            if (rangesToParse.Equals(""))
            {
                return;
            } 
            string[] ranges = rangesToParse.Split(new char[] {';'}, StringSplitOptions.RemoveEmptyEntries);
            for (int i = 0; i < ranges.Length; i++)
            {
                if (ranges[i].Contains("-"))
                {
                    string[] rangeParts = ranges[i].Split(new char[] {'-'}, StringSplitOptions.RemoveEmptyEntries);
                    if (rangeParts.Length != 2)
                    {
                        throw new HaploGrepException("Range could not be parsed to 2 elements "+ranges[i]);
                    }
                    int from = 0;
                    int to = 0;
                    try
                    {
                        from = (int)Convert.ToInt32(rangeParts[0].Trim());
                        to = (int)Convert.ToInt32(rangeParts[1].Trim());
                    }
                    catch (FormatException e1)
                    {
                        throw new HaploGrepException(rangeParts[0] + " " + rangeParts[1]);
                    } 
                    if ((to > 16570) || (from > 16570))
                    {
                        throw new HaploGrepException("Range specified was outside of valid rCRS range "+rangesToParse);
                    }
                    try
                    {
                        if (from > to)
                        {
                            addCustomRange(from, 16569);
                            addCustomRange(1, to);
                        }
                        else
                        {
                            addCustomRange(Convert.ToInt32(rangeParts[0].Trim()), Convert.ToInt32(rangeParts[1].Trim()));
                        }
                    }
                    catch (Exception e)
                    {
                        throw e;
                        //throw new HaploGrepException("Could not parse rangesToParse" +rangesToParse);
                    }
                }
                else
                {
                    try
                    {
                        addCustomRange(Convert.ToInt32(ranges[i].Trim()), Convert.ToInt32(ranges[i].Trim()));
                    }
                    catch (Exception e)
                    {
                        throw e;// new InvalidRangeException(rangesToParse);
                    }
                }
            }
        }
        public bool contains(Polymorphism polyToCheck)
        {
            for (int i = 0; i < this.starts.Count; i++)
            {
                if (((int)((int?)this.starts[i]) <= polyToCheck.Position) && ((int)((int?)this.ends[i]) >= polyToCheck.Position))
                {
                    return true;
                }
            } 
            return false;
        }
        public virtual void addCompleteRange()
        {
            this.starts.Add(Convert.ToInt32(1));
            this.ends.Add(Convert.ToInt32(16569));
        }
        public virtual void addCustomRange(int newRangeStart, int newRangeEnd)
        {
            this.starts.Add(Convert.ToInt32(newRangeStart));
            this.ends.Add(Convert.ToInt32(newRangeEnd));
        }
        public virtual string ToString()
        {
            string result = "";
            for (int i = 0; i < this.starts.Count; i++)
            {
                if (this.starts[i] == this.ends[i])
                {
                    result = result + this.starts[i] + " ; ";
                }
                else
                {
                    result = result + this.starts[i] + "-" + this.ends[i] + " ; ";
                }
            } 
            return result.Trim();
        }
    }
}