using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using HaploGrepSharp;
using System.Xml;

namespace HaploGrepSharp.TreeUtilities
{
    public class PhyloTreeNodev2 
    {
        public static int NodeIDCounter=0;
        public int NodeID;
        public static int NodesWithNoAccession = 0;
        public List<PhyloTreeNodev2> Children;
        public Haplogroup haplogroup;
        public PolymorphismCollection Mutations;
        public PhyloTreeNodev2(XmlNode node,PhyloTreeNodev2 parentNode=null)
        {
            this.NodeID=++NodeIDCounter;
            Children=new List<PhyloTreeNodev2>();
            Mutations=new PolymorphismCollection();
            //set the haplogroup
            haplogroup=new Haplogroup(node.Attributes.GetNamedItem("name").Value);
            var details=node.SelectSingleNode("details");
            haplogroup.AccessionId=details.Attributes.GetNamedItem("accessionNr").Value;
            if (String.IsNullOrEmpty(haplogroup.AccessionId))
            {
                NodesWithNoAccession++;
            }
            haplogroup.Reference=details.Attributes.GetNamedItem("reference").Value;
            //now copy polymorphism if needed
            if(parentNode!=null) {
                Mutations.AddRange(parentNode.Mutations);
            }
            //now update with the mutations here
            var polys =details.SelectNodes("poly");
            foreach(XmlNode p in polys)
            {
                if (p.InnerText.Contains("X"))
                    { //System.Console.WriteLine("Skipping: " + p.InnerText);
                    continue; }
                var currentPoly=new Polymorphism(p.InnerText);
                Mutations.Add(currentPoly);
            }
            //now make the children
            var children = node.SelectNodes("haplogroup");
            foreach (XmlNode currentElement in children)
            {
                Children.Add(new PhyloTreeNodev2(currentElement,this));
            }
        }
        public IEnumerable<PhyloTreeNodev2> GetAllChildren()
        {
            foreach (var node in Children)
            {
                yield return node;
                foreach (var subNode in node.GetAllChildren())
                {
                    yield return subNode;
                }
            }
        }
    }
       
}
