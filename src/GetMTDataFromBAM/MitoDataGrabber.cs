using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio;
using Bio.IO;
using Bio.IO.BAM;

namespace GetMTDataFromBAM
{
    /// <summary>
    /// This class contains utilities for taking a large BAM File and outputing the mitochondrial reads from that BAM File
    /// The matching regions it searchs are based on BLASTN to the entier genome, the python code to do this blast is given
    /// at the end of this source file.  The original file is called GetMTRegions.py
    /// </summary>
    public class MitoDataGrabber
    {
        public class MitoBlastMatchRegion
        {
            public readonly string Contig;
            public readonly int Start, End;
            public MitoBlastMatchRegion(string name, int start, int end)
            {
                //try
                //{
                //    Contig = "chr" + Int32.Parse(name).ToString();
                //}
                //catch
                //{
                   Contig = name;
                //}
                Start = start;
                End = end;
            }
        }
        #region GenomeAreasThatBlastToMTDNA
        public static readonly List<MitoBlastMatchRegion> regionsThatMatch = new List<MitoBlastMatchRegion>()
        {new MitoBlastMatchRegion("5",99390441,99381642),
    new MitoBlastMatchRegion("5",134264217,134258999),
    new MitoBlastMatchRegion("5",93903161,93906623),
    new MitoBlastMatchRegion("5",79948187,79945841),
    new MitoBlastMatchRegion("5",97747466,97745482),
    new MitoBlastMatchRegion("5",123096496,123097460),
    new MitoBlastMatchRegion("5",60057851,60057366),
    new MitoBlastMatchRegion("5",8622201,8622685),
    new MitoBlastMatchRegion("5",99390749,99390463),
    new MitoBlastMatchRegion("5",105889230,105889064),
    new MitoBlastMatchRegion("5",8620975,8621131),
    new MitoBlastMatchRegion("5",165957466,165957427),
    new MitoBlastMatchRegion("5",73071717,73071757),
    new MitoBlastMatchRegion("1",564464,570304),
    new MitoBlastMatchRegion("1",238109691,238104069),
    new MitoBlastMatchRegion("1",235705676,235700704),
    new MitoBlastMatchRegion("1",142792623,142793120),
    new MitoBlastMatchRegion("1",143344701,143345190),
    new MitoBlastMatchRegion("1",5910528,5910318),
    new MitoBlastMatchRegion("1",50482764,50483177),
    new MitoBlastMatchRegion("1",181391921,181392314),
    new MitoBlastMatchRegion("1",9634690,9634887),
    new MitoBlastMatchRegion("1",169443396,169443301),
    new MitoBlastMatchRegion("1",38077348,38077421),
    new MitoBlastMatchRegion("1",104163820,104163772),
    new MitoBlastMatchRegion("1",215673139,215673177),
    new MitoBlastMatchRegion("17",22022346,22031841),
    new MitoBlastMatchRegion("17",22018521,22020726),
    new MitoBlastMatchRegion("17",19501872,19503387),
    new MitoBlastMatchRegion("17",22021361,22022346),
    new MitoBlastMatchRegion("17",51183094,51183746),
    new MitoBlastMatchRegion("17",19505527,19506861),
    new MitoBlastMatchRegion("17",19504577,19505197),
    new MitoBlastMatchRegion("17",19503698,19504273),
    new MitoBlastMatchRegion("17",22020727,22021181),
    new MitoBlastMatchRegion("17",42075084,42075151),
    new MitoBlastMatchRegion("17",78591422,78591382),
    new MitoBlastMatchRegion("4",156387619,156372796),
    new MitoBlastMatchRegion("4",65476259,65472741),
    new MitoBlastMatchRegion("4",117218897,117221465),
    new MitoBlastMatchRegion("4",93622975,93623774),
    new MitoBlastMatchRegion("4",12642259,12641918),
    new MitoBlastMatchRegion("4",129002937,129002560),
    new MitoBlastMatchRegion("4",78929915,78929676),
    new MitoBlastMatchRegion("4",163342693,163342519),
    new MitoBlastMatchRegion("4",27732197,27732046),
    new MitoBlastMatchRegion("4",56194327,56194457),
    new MitoBlastMatchRegion("4",90652984,90653168),
    new MitoBlastMatchRegion("4",47774381,47774289),
    new MitoBlastMatchRegion("4",117217783,117218016),
    new MitoBlastMatchRegion("4",5406266,5406444),
    new MitoBlastMatchRegion("2",131029383,131040904),
    new MitoBlastMatchRegion("2",117778789,117784603),
    new MitoBlastMatchRegion("2",140981793,140977240),
    new MitoBlastMatchRegion("2",132141275,132137199),
    new MitoBlastMatchRegion("2",212641333,212639158),
    new MitoBlastMatchRegion("2",132143811,132141276),
    new MitoBlastMatchRegion("2",83042106,83046518),
    new MitoBlastMatchRegion("2",202079704,202077019),
    new MitoBlastMatchRegion("2",120970237,120973173),
    new MitoBlastMatchRegion("2",143849153,143847568),
    new MitoBlastMatchRegion("2",156121314,156119977),
    new MitoBlastMatchRegion("2",212641934,212643696),
    new MitoBlastMatchRegion("2",88124409,88124884),
    new MitoBlastMatchRegion("2",50816553,50815826),
    new MitoBlastMatchRegion("2",140975540,140974570),
    new MitoBlastMatchRegion("2",95565248,95564751),
    new MitoBlastMatchRegion("2",238432552,238431280),
    new MitoBlastMatchRegion("2",49457038,49456767),
    new MitoBlastMatchRegion("2",180604467,180604032),
    new MitoBlastMatchRegion("2",149639426,149639295),
    new MitoBlastMatchRegion("2",227586985,227587142),
    new MitoBlastMatchRegion("2",41012097,41012257),
    new MitoBlastMatchRegion("2",143850185,143849916),
    new MitoBlastMatchRegion("2",143846431,143846292),
    new MitoBlastMatchRegion("2",148022761,148022840),
    new MitoBlastMatchRegion("2",33992590,33992538),
    new MitoBlastMatchRegion("X",55210460,55204668),
    new MitoBlastMatchRegion("X",125605687,125606435),
    new MitoBlastMatchRegion("X",125606714,125607267),
    new MitoBlastMatchRegion("X",142517818,142519330),
    new MitoBlastMatchRegion("X",125606715,125606450),
    new MitoBlastMatchRegion("X",102064258,102064911),
    new MitoBlastMatchRegion("X",142519612,142519895),
    new MitoBlastMatchRegion("X",125864650,125864915),
    new MitoBlastMatchRegion("X",125865702,125865832),
    new MitoBlastMatchRegion("X",110660389,110660303),
    new MitoBlastMatchRegion("X",83659015,83658934),
    new MitoBlastMatchRegion("9",5092087,5098698),
    new MitoBlastMatchRegion("9",33656609,33659132),
    new MitoBlastMatchRegion("9",83180846,83178486),
    new MitoBlastMatchRegion("9",95301890,95301662),
    new MitoBlastMatchRegion("9",81425870,81425803),
    new MitoBlastMatchRegion("11",10531883,10529434),
    new MitoBlastMatchRegion("11",103280288,103272828),
    new MitoBlastMatchRegion("11",81263477,81268360),
    new MitoBlastMatchRegion("11",103281759,103280385),
    new MitoBlastMatchRegion("11",87524999,87524440),
    new MitoBlastMatchRegion("11",73221706,73221868),
    new MitoBlastMatchRegion("11",63954775,63954943),
    new MitoBlastMatchRegion("11",122874385,122874314),
    new MitoBlastMatchRegion("7",57257937,57265626),
    new MitoBlastMatchRegion("7",63572935,63568210),
    new MitoBlastMatchRegion("7",142373009,142375529),
    new MitoBlastMatchRegion("7",57237047,57241178),
    new MitoBlastMatchRegion("7",141504465,141501137),
    new MitoBlastMatchRegion("7",57253472,57255541),
    new MitoBlastMatchRegion("7",141505330,141504769),
    new MitoBlastMatchRegion("7",45291563,45291726),
    new MitoBlastMatchRegion("7",68798990,68798629),
    new MitoBlastMatchRegion("7",145694426,145694525),
    new MitoBlastMatchRegion("7",68201512,68201620),
    new MitoBlastMatchRegion("20",55935732,55932464),
    new MitoBlastMatchRegion("20",9149571,9149612),
    new MitoBlastMatchRegion("20",13147959,13148001),
    new MitoBlastMatchRegion("3",96337354,96336032),
    new MitoBlastMatchRegion("3",40293638,40295258),
    new MitoBlastMatchRegion("3",120440870,120441492),
    new MitoBlastMatchRegion("3",160665442,160665830),
    new MitoBlastMatchRegion("3",43271223,43270821),
    new MitoBlastMatchRegion("3",152637623,152637489),
    new MitoBlastMatchRegion("3",68708207,68708282),
    new MitoBlastMatchRegion("3",25508995,25509033),
    new MitoBlastMatchRegion("3",68719677,68719705),
    new MitoBlastMatchRegion("8",32868968,32872619),
    new MitoBlastMatchRegion("8",47742672,47739866),
    new MitoBlastMatchRegion("8",68499295,68496425),
    new MitoBlastMatchRegion("8",134767050,134768891),
    new MitoBlastMatchRegion("8",68496113,68493670),
    new MitoBlastMatchRegion("8",32872621,32873381),
    new MitoBlastMatchRegion("8",68500674,68499745),
    new MitoBlastMatchRegion("8",77113998,77114374),
    new MitoBlastMatchRegion("8",104095283,104096462),
    new MitoBlastMatchRegion("8",77557904,77557695),
    new MitoBlastMatchRegion("8",20408707,20408917),
    new MitoBlastMatchRegion("8",73897944,73898089),
    new MitoBlastMatchRegion("8",100508181,100508098),
    new MitoBlastMatchRegion("8",39928102,39928171),
    new MitoBlastMatchRegion("14",32954324,32953304),
    new MitoBlastMatchRegion("14",84637688,84639184),
    new MitoBlastMatchRegion("14",52054486,52054083),
    new MitoBlastMatchRegion("10",37891882,37889725),
    new MitoBlastMatchRegion("10",71353077,71350359),
    new MitoBlastMatchRegion("10",36724078,36721804),
    new MitoBlastMatchRegion("10",57358347,57359841),
    new MitoBlastMatchRegion("10",20036450,20037579),
    new MitoBlastMatchRegion("10",20035675,20036099),
    new MitoBlastMatchRegion("10",27162190,27162334),
    new MitoBlastMatchRegion("6",153986707,153990105),
    new MitoBlastMatchRegion("6",62284008,62284534),
    new MitoBlastMatchRegion("6",153990330,153990915),
    new MitoBlastMatchRegion("6",95157160,95156823),
    new MitoBlastMatchRegion("6",156869187,156868970),
    new MitoBlastMatchRegion("6",89265024,89265136),
    new MitoBlastMatchRegion("16",3421355,3419598),
    new MitoBlastMatchRegion("16",20733723,20732370),
    new MitoBlastMatchRegion("16",3422332,3421716),
    new MitoBlastMatchRegion("16",69392598,69392721),
    new MitoBlastMatchRegion("12",42093881,42092553),
    new MitoBlastMatchRegion("12",41757525,41757437),
    new MitoBlastMatchRegion("12",22158766,22158870),
    new MitoBlastMatchRegion("12",63167790,63167857),
    new MitoBlastMatchRegion("GL000211.1",79452,79949),
    new MitoBlastMatchRegion("GL000218.1",17026,17523),
    new MitoBlastMatchRegion("21",9735546,9736035),
    new MitoBlastMatchRegion("21",10492883,10493042),
    new MitoBlastMatchRegion("21",46796299,46796106),
    new MitoBlastMatchRegion("21",45895055,45894977),
    new MitoBlastMatchRegion("Y",13290173,13290670),
    new MitoBlastMatchRegion("Y",4212892,4212822),
    new MitoBlastMatchRegion("Y",8979505,8979570),
    new MitoBlastMatchRegion("13",110076727,110076472),
    new MitoBlastMatchRegion("13",56545890,56545768),
    new MitoBlastMatchRegion("13",41342484,41342558),
    new MitoBlastMatchRegion("18",45379808,45379617),
    new MitoBlastMatchRegion("18",59542118,59541786),
    new MitoBlastMatchRegion("18",2842352,2842198),
    new MitoBlastMatchRegion("22",36281719,36281765),
    new MitoBlastMatchRegion("15",67333249,67333291),
        new MitoBlastMatchRegion("MT",0,Int16.MaxValue)};
        #endregion

        public static IEnumerable<ISequence> OutputMitoReadsFromBamFile(string Filename)
        {
            //Very inefficient right now
            BAMParser bp = new BAMParser();
            foreach (MitoBlastMatchRegion region in regionsThatMatch)
            {
                    var map = bp.ParseRange(Filename, region.Contig, region.Start, region.End);
                    foreach (ISequence seq in map.QuerySequences.SelectMany(x => x.Sequences))
                    {
                        yield return seq;
                    }
            }
        }
        public static IEnumerable<ISequenceAlignmentParser> OutputAlignedMitoReadsFromBamFile(string Filename)
        {
            //Very inefficient right now
            BAMParser bp = new BAMParser();
            foreach (MitoBlastMatchRegion region in regionsThatMatch)
            {
                var map = bp.ParseRange(Filename, region.Contig, region.Start, region.End);
                foreach (ISequence seq in map.QuerySequences.SelectMany(x => x.Sequences))
                {
                    yield return seq as ISequenceAlignmentParser;
                }
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Filename"></param>
        /// <param name="fraction"></param>
        /// <param name="RemoveNames"></param>
        /// <returns></returns>
        public static IEnumerable<ISequence> OutputMitoReadsFromBamFileAlignedToCRSOnly(string Filename,double fraction,bool RemoveNames=true)
        {
            yield return null;
            //Random r = new Random();
            ////Very inefficient right now
            //BAMParser bp = new BAMParser();
            //foreach (MitoBlastMatchRegion region in regionsThatMatch)
            //{
            //    var map = bp.ParseRangeAsEnumerableSequences(Filename, "CRS", 0, Int16.MaxValue);
            //    foreach (var s in map)
            //    {
            //        if (RemoveNames)
            //        {
            //            s.QuerySequences[0].ID = "A";
            //        }
            //        if (fraction < 1)
            //        {
            //            if (r.NextDouble() < fraction)
            //                yield return s.QuerySequences[0];
            //        }
            //        else
            //        {
            //            yield return s.QuerySequences[0];
            //        }
            //    }
            //    //foreach (ISequence seq in map.QuerySequences.SelectMany(x => x.Sequences))
            //    //{
            //    //    yield return seq;
            //    //}
            //}
        }

    }
}

//import os
//from Bio import SeqIO
//refDirec=r"D:\MTData\MT_Data_Grabber\\"
//refFile=r"human_g1k_v37.fasta"
//queryFile=r"D:\MTData\MTNoN.fna"
//os.chdir(refDirec)
//infile=refFile
//cmd="makeblastdb -in="+infile+" -input_type=fasta -dbtype=nucl"
//print cmd
//os.system(cmd)
//print "LLL"
//cmd="blastn -query="+queryFile+" -outfmt=7 -out=BlastResult.txt -db="+infile
//print(cmd)
//os.system(cmd)

//outFile=open("BlastResult.txt")
//for line in outFile.readlines():
//    isinstance(line,str)
//    if line.startswith("#"):
//        continue
//    ls=line.split("\t")
//    chrm=ls[1]
//    if chrm=="MT":
//        continue
//    start=int(ls[8])
//    end=int(ls[9])
//    #print "chr"+chrm+":"+str(start)+"-"+str(end)
//    print 'new MitoBlastMatchRegion("'+chrm+'",'+str(start)+","+str(end)+"),"
// print 'new MitoBlastMatchRegion("MT",0,Int16.MaxValue)'