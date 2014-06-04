using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using HaploGrepSharp.NewSearchMethods;
using System.Diagnostics;

namespace HaploGrepSharp
{
	class Program
	{
		static string TestHSDFile = @"C:\Users\Delaney\SkyDrive\Software\MitochondrialPrograms\TestData\samplefile.hsd";

		static void Main (string[] args)
		{
			Tests.TestPhyloTreeWithGenbankData.RunTests ();
			//var i = new HaplotypeSearcher();
			//TreeUtilities.TreeLoader.OutputTreeAccessionNumbers("PhyloTreeAccessions.txt", CONSTANTS.TREE_XML_FILE);
			//HaploGrepIdentifier id = new HaploGrepIdentifier(TestHSDFile, "tmp.txt");

		}
	}
}
