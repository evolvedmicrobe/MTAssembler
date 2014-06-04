using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Xml;
using System.Linq;
using HaploGrepSharp;

namespace HaploGrepSharp.TreeUtilities
{
    public class TreeLoader
    {
        /// <summary>
        /// Load the phylotree XML file into a concrete object for use.
        /// </summary>
        /// <param name="xmlphyloTreeFile"></param>
        /// <returns></returns>
        public static PhyloTreeNodev2 LoadTree(string xmlphyloTreeFile)
        {
            //Arguments are the tree and the fluctation rates
            if (String.IsNullOrEmpty(xmlphyloTreeFile) || !File.Exists(xmlphyloTreeFile))
            {
                var fname = xmlphyloTreeFile ?? "No file set";
                throw new HaploGrepException("Could not load Tree XML File:" + fname);
            }
            try
            {
                var xmlDoc = new XmlDocument();
                xmlDoc.Load(xmlphyloTreeFile);
                XmlNode node = xmlDoc.SelectSingleNode("/phylotree").SelectSingleNode("haplogroup"); ;//.RootElement;
                return new PhyloTreeNodev2(node);
            }
            catch (Exception e)
            {
                throw e;
            }
        }

        public static void OutputTreeAccessionNumbers(string outFile, string xmlphylTreeFile)
        {
            var tree = LoadTree(xmlphylTreeFile);
            StreamWriter sw = new StreamWriter(outFile);
            foreach (var node in tree.GetAllChildren())
            {
                if (!String.IsNullOrEmpty(node.haplogroup.AccessionId))
                {
                    sw.WriteLine(node.haplogroup.AccessionId);
                }
            }
            sw.Close();
        }

    }
}
