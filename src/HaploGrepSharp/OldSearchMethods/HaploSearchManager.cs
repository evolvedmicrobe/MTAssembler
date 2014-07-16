using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Xml;
using System.Linq;

namespace HaploGrepSharp
{
    public sealed class HaploSearchManager
    {
        private List<Polymorphism> allPolysUsedinPhylotree = null;
        XmlDocument xmlDoc;
        public HaploSearchManager()
        {
            //Arguments are the tree and the fluctation rates
            this.allPolysUsedinPhylotree = new List<Polymorphism>();
            try {
                xmlDoc = new XmlDocument();
				var reader = new StringReader(PhyloTree15.TREE_STRING);
                xmlDoc.Load(reader);
                var skipped=SetPhylogeneticWeights();
                extractAllPolysFromPhylotree();
            }
            catch (Exception e) {
                throw e;
            }
        }
        private void extractAllPolysFromPhylotree()
        {
            XmlNodeList nameList = xmlDoc.SelectNodes("//poly");
            foreach (XmlNode a in nameList) {
                this.allPolysUsedinPhylotree.Add(new Polymorphism(a.InnerText));
            }
        }
        public static List<string> SetPhylogeneticWeights()
        {
            var skippedPositions = new List<string>();
            //StreamReader value = new StreamReader(fileName);
			var value = new StringReader(FlucRates15.FileAsString);
            string line = value.ReadLine();
            while (line != null) {
                string[] tokens = line.Split('\t');
                string polyString = tokens[0];
                double phyloGeneticWeight = Convert.ToDouble(tokens[1]);
                if (!polyString.Contains("X")) {
                    Polymorphism poly = new Polymorphism(polyString);
                    Polymorphism.changePhyloGeneticWeight(poly, phyloGeneticWeight);
                }
                else { 
                    //Console.WriteLine("Skipping: " + polyString);
                    skippedPositions.Add(polyString);
                }
                line = value.ReadLine();               
            }
            return skippedPositions;
        }
        public void changePoly(Haplogroup hg, Polymorphism polyOld, Polymorphism polyNew)
        {
            var e = getPolysOfHg(hg);
            foreach (XmlNode ce in e)
            {
                if (ce.InnerText.Equals(polyOld.ToString()))
                {
                    ce.InnerText = polyNew.ToString();
                    return;
                }
            }
            throw new HaploGrepException("Polymorphism does not exit");
        }
        public XmlDocument PhyloTree
        {
            get { return xmlDoc; }
        }
        public XmlNodeList getPolysOfHg(Haplogroup hg)
        {
            var titleNode = xmlDoc.SelectSingleNode("//haplogroup[@name=\"" + hg.ToString() + "\"]/details");
            return titleNode.SelectNodes("poly");
        }
        public List<Polymorphism> AllPolysUsedInPhylotree
        {
            get
            {
                return this.allPolysUsedinPhylotree;
            }
        }
    }
}