using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using System.Text.RegularExpressions;
namespace HaploGrepSharp
{
    public class Sample
    {
        public List<Polymorphism> sample = new List<Polymorphism>();
        public Sample(string sampleToParse, int callMethod)
        {
            var polyTokens = sampleToParse.Trim().Split(new char[]{' '},StringSplitOptions.RemoveEmptyEntries).ToList();
            this.sample = parseSample(polyTokens, callMethod);
        }
        public List<Polymorphism> Polymorphismn
        {
            get
            {
                return this.sample;
            }
        }
        public bool contains(Polymorphism polyToCheck)
        {
            return this.sample.Contains(polyToCheck);
        }
        public override string ToString()
        {
            string result = "";
            foreach (Polymorphism currentPoly in this.sample)
            {
                result = result + currentPoly + " ";
            }
            return result.Trim();
        }
        public void filter(List<Polymorphism> postivePolys)
        {
            List<Polymorphism> filteredSample = new List<Polymorphism>();
            foreach (Polymorphism currentPoly in this.sample)
            {
                if (postivePolys.Contains(currentPoly))
                {
                    filteredSample.Add(currentPoly);
                }
            }
            this.sample = filteredSample;
        }
        private Regex dotRegex=new Regex("\\d+\\.\\d\\w");
        private List<Polymorphism> parseSample(List<string> sample, int callMethod)
        {
            //List of e.g. G2a3	16183C	16189C	16193.1C	16223T	16278T	16362C	73G	260A	263G	309.1C	315.1C	489C	4833G	10400T	13563G	14569A	523d	524d		
            List<Polymorphism> filteredSample = new List<Polymorphism>();
            sample = sample.Select(x => x.Trim()).ToList();
            foreach (string currentPolyOrig in sample)
            {
                string currentPoly = currentPolyOrig;
                if ((!currentPoly.Contains("5899.1d!")) && (!currentPoly.Contains("65.1T(T)")))
                {
                    currentPoly = currentPoly.Replace("(", "");
                    currentPoly = currentPoly.Replace(")", "");
                    if (currentPoly.Contains("."))
                    {
                        if (callMethod == 0)
                        {
                            //309.1C
                            if (currentPoly.Contains(".1") && dotRegex.IsMatch(currentPoly))
                            {
                                string[] split = currentPoly.Split('.');
                                string position = split[0];//st1.nextToken();
                                string ins = split[1].Trim();// st1.nextToken().Trim();
                                int nextPos = 2;
                                //WTF? We are searching for another deletion at same location??
                                // Okay, really really really effing goofy here, explained by Sebastian below.
                                //If you enter it as 1235.1A, 1235.2A, 1235.3T in the hsd file, 
                                //Haplogrep transforms it to 1235.AAT.  including insertions in the current algorithm.
                                string samePos = position + ".";
                                foreach (string currentPoly2 in sample)
                                {
                                    if ((!currentPoly2.Equals(currentPoly)) && (currentPoly2.StartsWith(samePos)))
                                    {
                                        split = currentPoly2.Split('.');
                                        string position2 = split[0];
                                        string newIns = split[1].Trim();
                                        //Pattern p = Pattern.compile("\\d+");
                                        //Matcher m = p.matcher(ins);
                                        //m.find();
                                        //int tmp = Convert.ToInt32(ins.Substring(m.start(), m.end() - (m.start())));
                                        Match m = Regex.Match(newIns, "\\d+");
                                        int tmp = Convert.ToInt32(m.Value);
                                        if (nextPos != tmp)
                                        {
                                            throw new HaploGrepException(currentPoly2);
                                        }
                                        nextPos++;
                                        ins += newIns.Replace(m.Value, "");
                                    }
                                }
                                filteredSample.Add(new Polymorphism(position+"."+ins));
                            }
                                //Seems this will never be called as call method always 0
                            else if (Regex.IsMatch(currentPoly,"\\d+\\.1[a-zA-Z]{2,}")) {
                                Polymorphism newPoly = new Polymorphism(currentPoly);
                                filteredSample.Add(newPoly);
                            }
                        }
                        else {
                            Polymorphism newPoly = new Polymorphism(currentPoly);
                            filteredSample.Add(newPoly);
                        }
                    }
                    else if (currentPoly.Contains("-"))
                    {
                        string[] st1 = currentPoly.Split('-');
                        //StringTokenizer st1 = new StringTokenizer(currentPoly, "-");
                        string token = st1[0];// st1.nextToken();
                        int startPosition = (int)Convert.ToInt32(token);
                        token = st1[1];// st1.nextToken();
                        int endPosition = (int)Convert.ToInt32(token.Substring(0, token.Length - 1));
                        for (int i = startPosition; i <= endPosition; i++) {
                            filteredSample.Add(new Polymorphism(i, Mutations.DEL));
                            startPosition++;
                        }
                    }
                    else {
                        Polymorphism newPoly = new Polymorphism(currentPoly);
                        filteredSample.Add(newPoly);
                    }
                }
            }
            return filteredSample;
        }
    } 
}