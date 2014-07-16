﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio;
using Bio.Algorithms.MUMmer;
using Bio.Algorithms.Alignment;
using Bio.Extensions;

namespace HaploGrepSharp
{
    public static class ReferenceGenome
    {
        static NucMerQueryable nucmer;
        public static readonly Sequence ReferenceSequence;
        static ReferenceGenome()
        {
            ReferenceSequence = new Sequence(DnaAlphabet.Instance, rCRS, false);

            nucmer= new NucMerQueryable(ReferenceSequence);


        }
        public static IEnumerable<IEnumerable<DeltaAlignment>> GetDeltaAlignments(ISequence querySequence)
        {
            return nucmer.GetDeltaAlignments(querySequence);
            
        }
        /// <summary>
        /// Gets a reference sequence based on the 0 indexed values for position in the genom.e
        /// </summary>
        /// <param name="start">0 based index</param>
        /// <param name="end">0 based index, inclusive</param>
        /// <returns></returns>
        public static SectionOfReferenceGenome GetReferenceSequenceSection(int start, int end)
        {
            if (end > start)
            {
                var seq = ReferenceSequence.GetSubSequence(start, (end - start + 1));
                seq.ID = "Ref:" + start.ToString() + "-" + end.ToString();
                return new SectionOfReferenceGenome() { Start = start, End = end, Seq = seq as Sequence };
            }
            else
            {
                var seq1 = ReferenceSequence.GetSubSequence(start, ReferenceSequence.Count - start + 1);
                var seq2 = ReferenceSequence.GetSubSequence(0, end);
                List<byte> seqs = new List<byte>(seq1);
                seqs.AddRange(seq2);
                var seq = new Sequence(DnaAlphabet.Instance,seqs.ToArray());
                seq.ID = "Ref:" + start.ToString() + "-" + end.ToString();
                return new SectionOfReferenceGenome() {Start = start, End = end, Seq = seq};
            }

        }
        public const string DELTA_ALIGNMENT_METADATAKEY = "DeltaAlns";


		public static CompactSAMSequence AlignSequence(ISequence seq)
		{
            return null;
            //CompactSAMSequence css; 
            //if (seq is QualitativeSequence) {
            //    var qs = seq as QualitativeSequence;
            //    css = new CompactSAMSequence (seq.Alphabet,
            //        FastQFormatType.GATK_Recalibrated, 
            //        qs.ToArray (), 
            //        qs.GetPhredQualityScores ());
            //}
            //else
            //{
            //    css = new CompactSAMSequence (seq.Alphabet,
            //        FastQFormatType.GATK_Recalibrated, 
            //        seq.ToArray (), 
            //        Enumerable.Repeat ((byte)30, seq.Count).ToArray ());
            //}


		}

        public static void AssignContigToMTDNA(ISequence contig)
        {
            var delts = GetDeltaAlignments(contig).SelectMany(x => x).ToList();
            StringBuilder sb = new StringBuilder();
            int alnCounts = delts.Count;
            if (alnCounts > 0)
            {
                var alns = NucmerPairwiseAligner.ConvertDeltaToAlignment(delts);
                contig.Metadata[DELTA_ALIGNMENT_METADATAKEY] = alns;
                sb.Append("ALNS=" + alnCounts.ToString());                
                foreach (var d in alns)
                {
                    //if(d.IsReverseQueryDirection)
                    //{sb.Append("R;");}
                    //else{sb.Append("F;");}
                    var alnLength = d.FirstSequence.Count;
                    sb.Append(d.FirstOffset.ToString() + "-" + (d.FirstOffset+alnLength).ToString()+ ";" + d.SecondOffset.ToString() + "-" + (d.SecondOffset+alnLength).ToString() + ";" + d.Score.ToString());
                }
            }
            else {
                sb.Append("No Alignments");
            }
            contig.ID = contig.ID + " " + sb.ToString();
        }
        /// <summary>
        /// Converts a 0 based position on the reference sequence without the N to
        /// a one base position on the sequence with it.
        /// The rCRS has a "N" at position 3106 (0 based) or 3107 (1 based)
        /// </summary>
        /// <param name="pos">0 based position</param>
        /// <returns>adjusted value in one based index</returns>
        public static int ConvertTorCRSPosition(int pos)
        {
            return pos >= 3106 ? pos + 2 : pos + 1;
        }

		/// <summary>
		/// Gets the reference base at the 1-based index position.
		/// </summary>
		/// <param name="position"></param>
		/// <returns></returns>
		public static char GetReferenceBaseAt_rCRSPosition(int position)
		{
			int n_pos = position - 1;
			if (n_pos == 3106) { return 'N'; }
			n_pos = n_pos > 3106 ? n_pos-1 : n_pos;
			return rCRS[n_pos];

		}

        /// <summary>
        /// The reference cambridge sequence with the "N" removed;
        /// </summary>
        public readonly static string rCRS = @"gatcacaggtctatcaccctattaaccactcacgggagctctccatgcatttggtatttt
cgtctggggggtatgcacgcgatagcattgcgagacgctggagccggagcaccctatgtc
gcagtatctgtctttgattcctgcctcatcctattatttatcgcacctacgttcaatatt
acaggcgaacatacttactaaagtgtgttaattaattaatgcttgtaggacataataata
acaattgaatgtctgcacagccactttccacacagacatcataacaaaaaatttccacca
aaccccccctcccccgcttctggccacagcacttaaacacatctctgccaaaccccaaaa
acaaagaaccctaacaccagcctaaccagatttcaaattttatcttttggcggtatgcac
ttttaacagtcaccccccaactaacacattattttcccctcccactcccatactactaat
ctcatcaatacaacccccgcccatcctacccagcacacacacaccgctgctaaccccata
ccccgaaccaaccaaaccccaaagacaccccccacagtttatgtagcttacctcctcaaa
gcaatacactgaaaatgtttagacgggctcacatcaccccataaacaaataggtttggtc
ctagcctttctattagctcttagtaagattacacatgcaagcatccccgttccagtgagt
tcaccctctaaatcaccacgatcaaaaggaacaagcatcaagcacgcagcaatgcagctc
aaaacgcttagcctagccacacccccacgggaaacagcagtgattaacctttagcaataa
acgaaagtttaactaagctatactaaccccagggttggtcaatttcgtgccagccaccgc
ggtcacacgattaacccaagtcaatagaagccggcgtaaagagtgttttagatcaccccc
tccccaataaagctaaaactcacctgagttgtaaaaaactccagttgacacaaaatagac
tacgaaagtggctttaacatatctgaacacacaatagctaagacccaaactgggattaga
taccccactatgcttagccctaaacctcaacagttaaatcaacaaaactgctcgccagaa
cactacgagccacagcttaaaactcaaaggacctggcggtgcttcatatccctctagagg
agcctgttctgtaatcgataaaccccgatcaacctcaccacctcttgctcagcctatata
ccgccatcttcagcaaaccctgatgaaggctacaaagtaagcgcaagtacccacgtaaag
acgttaggtcaaggtgtagcccatgaggtggcaagaaatgggctacattttctaccccag
aaaactacgatagcccttatgaaacttaagggtcgaaggtggatttagcagtaaactaag
agtagagtgcttagttgaacagggccctgaagcgcgtacacaccgcccgtcaccctcctc
aagtatacttcaaaggacatttaactaaaacccctacgcatttatatagaggagacaagt
cgtaacatggtaagtgtactggaaagtgcacttggacgaaccagagtgtagcttaacaca
aagcacccaacttacacttaggagatttcaacttaacttgaccgctctgagctaaaccta
gccccaaacccactccaccttactaccagacaaccttagccaaaccatttacccaaataa
agtataggcgatagaaattgaaacctggcgcaatagatatagtaccgcaagggaaagatg
aaaaattataaccaagcataatatagcaaggactaacccctataccttctgcataatgaa
ttaactagaaataactttgcaaggagagccaaagctaagacccccgaaaccagacgagct
acctaagaacagctaaaagagcacacccgtctatgtagcaaaatagtgggaagatttata
ggtagaggcgacaaacctaccgagcctggtgatagctggttgtccaagatagaatcttag
ttcaactttaaatttgcccacagaaccctctaaatccccttgtaaatttaactgttagtc
caaagaggaacagctctttggacactaggaaaaaaccttgtagagagagtaaaaaattta
acacccatagtaggcctaaaagcagccaccaattaagaaagcgttcaagctcaacaccca
ctacctaaaaaatcccaaacatataactgaactcctcacacccaattggaccaatctatc
accctatagaagaactaatgttagtataagtaacatgaaaacattctcctccgcataagc
ctgcgtcagattaaaacactgaactgacaattaacagcccaatatctacaatcaaccaac
aagtcattattaccctcactgtcaacccaacacaggcatgctcataaggaaaggttaaaa
aaagtaaaaggaactcggcaaatcttaccccgcctgtttaccaaaaacatcacctctagc
atcaccagtattagaggcaccgcctgcccagtgacacatgtttaacggccgcggtaccct
aaccgtgcaaaggtagcataatcacttgttccttaaatagggacctgtatgaatggctcc
acgagggttcagctgtctcttacttttaaccagtgaaattgacctgcccgtgaagaggcg
ggcataacacagcaagacgagaagaccctatggagctttaatttattaatgcaaacagta
cctaacaaacccacaggtcctaaactaccaaacctgcattaaaaatttcggttggggcga
cctcggagcagaacccaacctccgagcagtacatgctaagacttcaccagtcaaagcgaa
ctactatactcaattgatccaataacttgaccaacggaacaagttaccctagggataaca
gcgcaatcctattctagagtccatatcaacaatagggtttacgacctcgatgttggatca
ggacatcccgatggtgcagccgctattaaaggttcgtttgttcaacgattaaagtcctac
gtgatctgagttcagaccggagtaatccaggtcggtttctatctacttcaaattcctcc
ctgtacgaaaggacaagagaaataaggcctacttcacaaagcgccttcccccgtaaatga
tatcatctcaacttagtattatacccacacccacccaagaacagggtttgttaagatggc
agagcccggtaatcgcataaaacttaaaactttacagtcagaggttcaattcctcttctt
aacaacatacccatggccaacctcctactcctcattgtacccattctaatcgcaatggca
ttcctaatgcttaccgaacgaaaaattctaggctatatacaactacgcaaaggccccaac
gttgtaggcccctacgggctactacaacccttcgctgacgccataaaactcttcaccaaa
gagcccctaaaacccgccacatctaccatcaccctctacatcaccgccccgaccttagct
ctcaccatcgctcttctactatgaacccccctccccatacccaaccccctggtcaacctc
aacctaggcctcctatttattctagccacctctagcctagccgtttactcaatcctctga
tcagggtgagcatcaaactcaaactacgccctgatcggcgcactgcgagcagtagcccaa
acaatctcatatgaagtcaccctagccatcattctactatcaacattactaataagtggc
tcctttaacctctccacccttatcacaacacaagaacacctctgattactcctgccatca
tgacccttggccataatatgatttatctccacactagcagagaccaaccgaacccccttc
gaccttgccgaaggggagtccgaactagtctcaggcttcaacatcgaatacgccgcaggc
cccttcgccctattcttcatagccgaatacacaaacattattataataaacaccctcacc
actacaatcttcctaggaacaacatatgacgcactctcccctgaactctacacaacatat
tttgtcaccaagaccctacttctaacctccctgttcttatgaattcgaacagcatacccc
cgattccgctacgaccaactcatacacctcctatgaaaaaacttcctaccactcacccta
gcattacttatatgatatgtctccatacccattacaatctccagcattccccctcaaacc
taagaaatatgtctgataaaagagttactttgatagagtaaataataggagcttaaaccc
ccttatttctaggactatgagaatcgaacccatccctgagaatccaaaattctccgtgcc
acctatcacaccccatcctaaagtaaggtcagctaaataagctatcgggcccataccccg
aaaatgttggttatacccttcccgtactaattaatcccctggcccaacccgtcatctact
ctaccatctttgcaggcacactcatcacagcgctaagctcgcactgattttttacctgag
taggcctagaaataaacatgctagcttttattccagttctaaccaaaaaaataaaccctc
gttccacagaagctgccatcaagtatttcctcacgcaagcaaccgcatccataatccttc
taatagctatcctcttcaacaatatactctccggacaatgaaccataaccaatactacca
atcaatactcatcattaataatcataatagctatagcaataaaactaggaatagccccct
ttcacttctgagtcccagaggttacccaaggcacccctctgacatccggcctgcttcttc
tcacatgacaaaaactagcccccatctcaatcatataccaaatctctccctcactaaacg
taagccttctcctcactctctcaatcttatccatcatagcaggcagttgaggtggattaa
accaaacccagctacgcaaaatcttagcatactcctcaattacccacataggatgaataa
tagcagttctaccgtacaaccctaacataaccattcttaatttaactatttatattatcc
taactactaccgcattcctactactcaacttaaactccagcaccacgaccctactactat
ctcgcacctgaaacaagctaacatgactaacacccttaattccatccaccctcctctccc
taggaggcctgcccccgctaaccggctttttgcccaaatgggccattatcgaagaattca
caaaaaacaatagcctcatcatccccaccatcatagccaccatcaccctccttaacctct
acttctacctacgcctaatctactccacctcaatcacactactccccatatctaacaacg
taaaaataaaatgacagtttgaacatacaaaacccaccccattcctccccacactcatcg
cccttaccacgctactcctacctatctccccttttatactaataatcttatagaaattta
ggttaaatacagaccaagagccttcaaagccctcagtaagttgcaatacttaatttctgt
aacagctaaggactgcaaaaccccactctgcatcaactgaacgcaaatcagccactttaa
ttaagctaagcccttactagaccaatgggacttaaacccacaaacacttagttaacagct
aagcaccctaatcaactggcttcaatctacttctcccgccgccgggaaaaaaggcgggag
aagccccggcaggtttgaagctgcttcttcgaatttgcaattcaatatgaaaatcacctc
ggagctggtaaaaagaggcctaacccctgtctttagatttacagtccaatgcttcactca
gccattttacctcacccccactgatgttcgccgaccgttgactattctctacaaaccaca
aagacattggaacactatacctattattcggcgcatgagctggagtcctaggcacagctc
taagcctccttattcgagccgagctgggccagccaggcaaccttctaggtaacgaccaca
tctacaacgttatcgtcacagcccatgcatttgtaataatcttcttcatagtaataccca
tcataatcggaggctttggcaactgactagttcccctaataatcggtgcccccgatatgg
cgtttccccgcataaacaacataagcttctgactcttacctccctctctcctactcctgc
tcgcatctgctatagtggaggccggagcaggaacaggttgaacagtctaccctcccttag
cagggaactactcccaccctggagcctccgtagacctaaccatcttctccttacacctag
caggtgtctcctctatcttaggggccatcaatttcatcacaacaattatcaatataaaac
cccctgccataacccaataccaaacgcccctcttcgtctgatccgtcctaatcacagcag
tcctacttctcctatctctcccagtcctagctgctggcatcactatactactaacagacc
gcaacctcaacaccaccttcttcgaccccgccggaggaggagaccccattctataccaac
acctattctgatttttcggtcaccctgaagtttatattcttatcctaccaggcttcggaa
taatctcccatattgtaacttactactccggaaaaaaagaaccatttggatacataggta
tggtctgagctatgatatcaattggcttcctagggtttatcgtgtgagcacaccatatat
ttacagtaggaatagacgtagacacacgagcatatttcacctccgctaccataatcatcg
ctatccccaccggcgtcaaagtatttagctgactcgccacactccacggaagcaatatga
aatgatctgctgcagtgctctgagccctaggattcatctttcttttcaccgtaggtggcc
tgactggcattgtattagcaaactcatcactagacatcgtactacacgacacgtactacg
ttgtagcccacttccactatgtcctatcaataggagctgtatttgccatcataggaggct
tcattcactgatttcccctattctcaggctacaccctagaccaaacctacgccaaaatcc
atttcactatcatattcatcggcgtaaatctaactttcttcccacaacactttctcggcc
tatccggaatgccccgacgttactcggactaccccgatgcatacaccacatgaaacatcc
tatcatctgtaggctcattcatttctctaacagcagtaatattaataattttcatgattt
gagaagccttcgcttcgaagcgaaaagtcctaatagtagaagaaccctccataaacctgg
agtgactatatggatgccccccaccctaccacacattcgaagaacccgtatacataaaat
ctagacaaaaaaggaaggaatcgaaccccccaaagctggtttcaagccaaccccatggcc
tccatgactttttcaaaaaggtattagaaaaaccatttcataactttgtcaaagttaaat
tataggctaaatcctatatatcttaatggcacatgcagcgcaagtaggtctacaagacgc
tacttcccctatcatagaagagcttatcacctttcatgatcacgccctcataatcatttt
ccttatctgcttcctagtcctgtatgcccttttcctaacactcacaacaaaactaactaa
tactaacatctcagacgctcaggaaatagaaaccgtctgaactatcctgcccgccatcat
cctagtcctcatcgccctcccatccctacgcatcctttacataacagacgaggtcaacga
tccctcccttaccatcaaatcaattggccaccaatggtactgaacctacgagtacaccga
ctacggcggactaatcttcaactcctacatacttcccccattattcctagaaccaggcga
cctgcgactccttgacgttgacaatcgagtagtactcccgattgaagcccccattcgtat
aataattacatcacaagacgtcttgcactcatgagctgtccccacattaggcttaaaaac
agatgcaattcccggacgtctaaaccaaaccactttcaccgctacacgaccgggggtata
ctacggtcaatgctctgaaatctgtggagcaaaccacagtttcatgcccatcgtcctaga
attaattcccctaaaaatctttgaaatagggcccgtatttaccctatagcaccccctcta
ccccctctagagcccactgtaaagctaacttagcattaaccttttaagttaaagattaag
agaaccaacacctctttacagtgaaatgccccaactaaatactaccgtatggcccaccat
aattacccccatactccttacactattcctcatcacccaactaaaaatattaaacacaaa
ctaccacctacctccctcaccaaagcccataaaaataaaaaattataacaaaccctgaga
accaaaatgaacgaaaatctgttcgcttcattcattgcccccacaatcctaggcctaccc
gccgcagtactgatcattctatttccccctctattgatccccacctccaaatatctcatc
aacaaccgactaatcaccacccaacaatgactaatcaaactaacctcaaaacaaatgata
accatacacaacactaaaggacgaacctgatctcttatactagtatccttaatcattttt
attgccacaactaacctcctcggactcctgcctcactcatttacaccaaccacccaacta
tctataaacctagccatggccatccccttatgagcgggcacagtgattataggctttcgc
tctaagattaaaaatgccctagcccacttcttaccacaaggcacacctacaccccttatc
cccatactagttattatcgaaaccatcagcctactcattcaaccaatagccctggccgta
cgcctaaccgctaacattactgcaggccacctactcatgcacctaattggaagcgccacc
ctagcaatatcaaccattaaccttccctctacacttatcatcttcacaattctaattcta
ctgactatcctagaaatcgctgtcgccttaatccaagcctacgttttcacacttctagta
agcctctacctgcacgacaacacataatgacccaccaatcacatgcctatcatatagtaa
aacccagcccatgacccctaacaggggccctctcagccctcctaatgacctccggcctag
ccatgtgatttcacttccactccataacgctcctcatactaggcctactaaccaacacac
taaccatataccaatgatggcgcgatgtaacacgagaaagcacataccaaggccaccaca
caccacctgtccaaaaaggccttcgatacgggataatcctatttattacctcagaagttt
ttttcttcgcaggatttttctgagccttttaccactccagcctagcccctaccccccaat
taggagggcactggcccccaacaggcatcaccccgctaaatcccctagaagtcccactcc
taaacacatccgtattactcgcatcaggagtatcaatcacctgagctcaccatagtctaa
tagaaaacaaccgaaaccaaataattcaagcactgcttattacaattttactgggtctct
attttaccctcctacaagcctcagagtacttcgagtctcccttcaccatttccgacggca
tctacggctcaacattttttgtagccacaggcttccacggacttcacgtcattattggct
caactttcctcactatctgcttcatccgccaactaatatttcactttacatccaaacatc
actttggcttcgaagccgccgcctgatactggcattttgtagatgtggtttgactatttc
tgtatgtctccatctattgatgagggtcttactcttttagtataaatagtaccgttaact
tccaattaactagttttgacaacattcaaaaaagagtaataaacttcgccttaattttaa
taatcaacaccctcctagccttactactaataattattacattttgactaccacaactca
acggctacatagaaaaatccaccccttacgagtgcggcttcgaccctatatcccccgccc
gcgtccctttctccataaaattcttcttagtagctattaccttcttattatttgatctag
aaattgccctccttttacccctaccatgagccctacaaacaactaacctgccactaatag
ttatgtcatccctcttattaatcatcatcctagccctaagtctggcctatgagtgactac
aaaaaggattagactgaaccgaattggtatatagtttaaacaaaacgaatgatttcgact
cattaaattatgataatcatatttaccaaatgcccctcatttacataaatattatactag
catttaccatctcacttctaggaatactagtatatcgctcacacctcatatcctccctac
tatgcctagaaggaataatactatcgctgttcattatagctactctcataaccctcaaca
cccactccctcttagccaatattgtgcctattgccatactagtctttgccgcctgcgaag
cagcggtgggcctagccctactagtctcaatctccaacacatatggcctagactacgtac
ataacctaaacctactccaatgctaaaactaatcgtcccaacaattatattactaccact
gacatgactttccaaaaaacacataatttgaatcaacacaaccacccacagcctaattat
tagcatcatccctctactattttttaaccaaatcaacaacaacctatttagctgttcccc
aaccttttcctccgaccccctaacaacccccctcctaatactaactacctgactcctacc
cctcacaatcatggcaagccaacgccacttatccagtgaaccactatcacgaaaaaaact
ctacctctctatactaatctccctacaaatctccttaattataacattcacagccacaga
actaatcatattttatatcttcttcgaaaccacacttatccccaccttggctatcatcac
ccgatgaggcaaccagccagaacgcctgaacgcaggcacatacttcctattctacaccct
agtaggctcccttcccctactcatcgcactaatttacactcacaacaccctaggctcact
aaacattctactactcactctcactgcccaagaactatcaaactcctgagccaacaactt
aatatgactagcttacacaatagcttttatagtaaagatacctctttacggactccactt
atgactccctaaagcccatgtcgaagcccccatcgctgggtcaatagtacttgccgcagt
actcttaaaactaggcggctatggtataatacgcctcacactcattctcaaccccctgac
aaaacacatagcctaccccttccttgtactatccctatgaggcataattataacaagctc
catctgcctacgacaaacagacctaaaatcgctcattgcatactcttcaatcagccacat
agccctcgtagtaacagccattctcatccaaaccccctgaagcttcaccggcgcagtcat
tctcataatcgcccacgggcttacatcctcattactattctgcctagcaaactcaaacta
cgaacgcactcacagtcgcatcataatcctctctcaaggacttcaaactctactcccact
aatagctttttgatgacttctagcaagcctcgctaacctcgccttaccccccactattaa
cctactgggagaactctctgtgctagtaaccacgttctcctgatcaaatatcactctcct
acttacaggactcaacatactagtcacagccctatactccctctacatatttaccacaac
acaatggggctcactcacccaccacattaacaacataaaaccctcattcacacgagaaaa
caccctcatgttcatacacctatcccccattctcctcctatccctcaaccccgacatcat
taccgggttttcctcttgtaaatatagtttaaccaaaacatcagattgtgaatctgacaa
cagaggcttacgaccccttatttaccgagaaagctcacaagaactgctaactcatgcccc
catgtctaacaacatggctttctcaacttttaaaggataacagctatccattggtcttag
gccccaaaaattttggtgcaactccaaataaaagtaataaccatgcacactactataacc
accctaaccctgacttccctaattccccccatccttaccaccctcgttaaccctaacaaa
aaaaactcatacccccattatgtaaaatccattgtcgcatccacctttattatcagtctc
ttccccacaacaatattcatgtgcctagaccaagaagttattatctcgaactgacactga
gccacaacccaaacaacccagctctccctaagcttcaaactagactacttctccataata
ttcatccctgtagcattgttcgttacatggtccatcatagaattctcactgtgatatata
aactcagacccaaacattaatcagttcttcaaatatctactcatcttcctaattaccata
ctaatcttagttaccgctaacaacctattccaactgttcatcggctgagagggcgtagga
attatatccttcttgctcatcagttgatgatacgcccgagcagatgccaacacagcagcc
attcaagcaatcctatacaaccgtatcggcgatatcggtttcatcctcgccttagcatga
tttatcctacactccaactcatgagacccacaacaaatagcccttctaaacgctaatcca
agcctcaccccactactaggcctcctcctagcagcagcaggcaaatcagcccaattaggt
ctccacccctgactcccctcagccatagaaggccccaccccagtctcagccctactccac
tcaagcactatagttgtagcaggaatcttcttactcatccgcttccaccccctagcagaa
aatagcccactaatccaaactctaacactatgcttaggcgctatcaccactctgttcgca
gcagtctgcgcccttacacaaaatgacatcaaaaaaatcgtagccttctccacttcaagt
caactaggactcataatagttacaatcggcatcaaccaaccacacctagcattcctgcac
atctgtacccacgccttcttcaaagccatactatttatgtgctccgggtccatcatccac
aaccttaacaatgaacaagatattcgaaaaataggaggactactcaaaaccatacctctc
acttcaacctccctcaccattggcagcctagcattagcaggaatacctttcctcacaggt
ttctactccaaagaccacatcatcgaaaccgcaaacatatcatacacaaacgcctgagcc
ctatctattactctcatcgctacctccctgacaagcgcctatagcactcgaataattctt
ctcaccctaacaggtcaacctcgcttccccacccttactaacattaacgaaaataacccc
accctactaaaccccattaaacgcctggcagccggaagcctattcgcaggatttctcatt
actaacaacatttcccccgcatcccccttccaaacaacaatccccctctacctaaaactc
acagccctcgctgtcactttcctaggacttctaacagccctagacctcaactacctaacc
aacaaacttaaaataaaatccccactatgcacattttatttctccaacatactcggattc
taccctagcatcacacaccgcacaatcccctatctaggccttcttacgagccaaaacctg
cccctactcctcctagacctaacctgactagaaaagctattacctaaaacaatttcacag
caccaaatctccacctccatcatcacctcaacccaaaaaggcataattaaactttacttc
ctctctttcttcttcccactcatcctaaccctactcctaatcacataacctattcccccg
agcaatctcaattacaatatatacaccaacaaacaatgttcaaccagtaactactactaa
tcaacgcccataatcatacaaagcccccgcaccaataggatcctcccgaatcaaccctga
cccctctccttcataaattattcagcttcctacactattaaagtttaccacaaccaccac
cccatcatactctttcacccacagcaccaatcctacctccatcgctaaccccactaaaac
actcaccaagacctcaacccctgacccccatgcctcaggatactcctcaatagccatcgc
tgtagtatatccaaagacaaccatcattccccctaaataaattaaaaaaactattaaacc
catataacctcccccaaaattcagaataataacacacccgaccacaccgctaacaatcaa
tactaaacccccataaataggagaaggcttagaagaaaaccccacaaaccccattactaa
acccacactcaacagaaacaaagcatacatcattattctcgcacggactacaaccacgac
caatgatatgaaaaaccatcgttgtatttcaactacaagaacaccaatgaccccaatacg
caaaactaaccccctaataaaattaattaaccactcattcatcgacctccccaccccatc
caacatctccgcatgatgaaacttcggctcactccttggcgcctgcctgatcctccaaat
caccacaggactattcctagccatgcactactcaccagacgcctcaaccgccttttcatc
aatcgcccacatcactcgagacgtaaattatggctgaatcatccgctaccttcacgccaa
tggcgcctcaatattctttatctgcctcttcctacacatcgggcgaggcctatattacgg
atcatttctctactcagaaacctgaaacatcggcattatcctcctgcttgcaactatagc
aacagccttcataggctatgtcctcccgtgaggccaaatatcattctgaggggccacagt
aattacaaacttactatccgccatcccatacattgggacagacctagttcaatgaatctg
aggaggctactcagtagacagtcccaccctcacacgattctttacctttcacttcatctt
gcccttcattattgcagccctagcaacactccacctcctattcttgcacgaaacgggatc
aaacaaccccctaggaatcacctcccattccgataaaatcaccttccacccttactacac
aatcaaagacgccctcggcttacttctcttccttctctccttaatgacattaacactatt
ctcaccagacctcctaggcgacccagacaattataccctagccaaccccttaaacacccc
tccccacatcaagcccgaatgatatttcctattcgcctacacaattctccgatccgtccc
taacaaactaggaggcgtccttgccctattactatccatcctcatcctagcaataatccc
catcctccatatatccaaacaacaaagcataatatttcgcccactaagccaatcacttta
ttgactcctagccgcagacctcctcattctaacctgaatcggaggacaaccagtaagcta
cccttttaccatcattggacaagtagcatccgtactatacttcacaacaatcctaatcct
aataccaactatctccctaattgaaaacaaaatactcaaatgggcctgtccttgtagtat
aaactaatacaccagtcttgtaaaccggagatgaaaacctttttccaaggacaaatcaga
gaaaaagtctttaactccaccattagcacccaaagctaagattctaatttaaactattct
ctgttctttcatggggaagcagatttgggtaccacccaagtattgactcacccatcaaca
accgctatgtatttcgtacattactgccagccaccatgaatattgtacggtaccataaat
acttgaccacctgtagtacataaaaacccaatccacatcaaaaccccctccccatgctta
caagcaagtacagcaatcaaccctcaactatcacacatcaactgcaactccaaagccacc
cctcacccactaggataccaacaaacctacccacccttaacagtacatagtacataaagc
catttaccgtacatagcacattacagtcaaatcccttctcgtccccatggatgacccccc
tcagataggggtcccttgaccaccatcctccgtgaaatcaatatcccgcacaagagtgct
actctcctcgctccgggcccataacacttgggggtagctaaagtgaactgtatccgacat
ctggttcctacttcagggtcataaagcctaaatagcccacacgttccccttaaataagac
atcacgatg".Replace("\n","").Replace("\r","").ToUpper();


    }
}
