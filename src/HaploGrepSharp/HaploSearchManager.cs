using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Xml;
using System.Linq;

namespace HaploGrep
{
    public sealed class HaploSearchManager
    {
        private List<Polymorphism> allPolysUsedinPhylotree = null;
        private string phyloString;
        XmlDocument xmlDoc;
        public HaploSearchManager(string phylotree, string weights)
        {
            //Arguments are the tree and the fluctation rates
            this.allPolysUsedinPhylotree = new List<Polymorphism>();
            this.phyloString = phylotree;
			if (String.IsNullOrEmpty(phylotree) || !File.Exists (phylotree)) {
				throw new HaploGrepException ("Could not find file: " + phylotree);
			}
			if (String.IsNullOrEmpty (weights) || !File.Exists (weights)) {
				throw new HaploGrepException ("Could not find file: " + weights);
			}
            try
            {
                xmlDoc = new XmlDocument();
                xmlDoc.Load(phylotree);
                LoadHaplotypes();
                //Files phyloFile = this.GetType().ClassLoader.getResourceAsStream(phylotree);
                //InputStream flucRates = this.GetType().ClassLoader.getResourceAsStream(weights);
                //this.phyloTree = builder.build(phyloFile);
				SetPolygeneticWeights(new StreamReader(weights));
                extractAllPolysFromPhylotree();
            }
            catch (Exception e)
            {
				throw new HaploGrepException("Could not load HaploSearchManager.", e);
            }
        }
        void LoadHaplotypes()
        {
            XmlNodeList xnd = xmlDoc.SelectNodes("//haplogroup");
            List<Haplogroup> hg = new List<Haplogroup>();
            foreach (XmlNode n in xnd)
            {
                var details = n.SelectSingleNode("/details");
            }
        }
        public string PhyloString
        {
			get {
                return this.phyloString;
            }
            set {
                this.phyloString = value;
            }
        }
        private void extractAllPolysFromPhylotree()
        {
            XmlNodeList nameList = xmlDoc.SelectNodes("//poly");
            foreach (XmlNode a in nameList)
            {
                this.allPolysUsedinPhylotree.Add(new Polymorphism(a.InnerText));
            }
        }
		private void SetPolygeneticWeights(StreamReader inFile)
		{
			string line = inFile.ReadLine();
			while (line != null)
			{
				string[] tokens = line.Split('\t');
				string polyString = tokens[0];
				double phyloGeneticWeight = Convert.ToDouble(tokens[1]);
				if (!polyString.Contains("X"))
				{
					Polymorphism poly = new Polymorphism(polyString);
					Polymorphism.changePhyloGeneticWeight(poly, this.phyloString, phyloGeneticWeight);
				}
				else { Console.WriteLine("Skipping: "+polyString); }
				line = inFile.ReadLine();                    
			}
		}

		//TODO: This needs to be a method
		//private StreamReader PolygeneticWeights {


		public void changePoly(Haplogroup hg, Polymorphism polyOld, Polymorphism polyNew) {
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
		public XmlDocument PhyloTree {
            get { return xmlDoc; }
        }
		public XmlNodeList getPolysOfHg(Haplogroup hg) {
            var titleNode = xmlDoc.SelectSingleNode("//haplogroup[@name=\"" + hg.ToString() + "\"]/details");
            return titleNode.SelectNodes("poly");
        }
		public List<Polymorphism> AllPolysUsedInPhylotree {
			get {
                return this.allPolysUsedinPhylotree;
            }
        }
    }
}