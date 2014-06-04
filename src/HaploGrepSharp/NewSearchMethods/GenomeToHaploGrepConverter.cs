using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio;

namespace HaploGrepSharp.NewSearchMethods
{
    /// <summary>
    /// Takes a genome sequence of a full mitochondria, aligns it and produces a list of polymorphisms in haplogrep format
    /// </summary>
    public class GenomeToHaploGrepConverter
    {
        /// <summary>
        /// Take a mtDNA genome and align it to the reference genome, producing
        /// all polymorphisms found
        /// </summary>
        /// <param name="toConvert"></param>
        /// <returns></returns>
        public static List<Polymorphism> FindPolymorphisms(Sequence toConvert)
        {
            var delts = ReferenceGenome.GetDeltaAlignments(toConvert).SelectMany(x => x).ToList();
            if (delts.Count != 1) {
                throw new HaploGrepException("Genome assembly had no or multiple delta alignments with the reference, whereas only 1 is expected!");
            }           
            var delt = delts[0];
            var aln = delt.ConvertDeltaToSequences();
            var refseq = aln.FirstSequence as Sequence;
            var qseq = aln.SecondSequence as Sequence;

            int position = (int)delt.FirstSequenceStart;
            List<Polymorphism> res = new List<Polymorphism>();
            for (int i = 0; i < refseq.Count; i++) {
                if (refseq[i] != qseq[i])
                {
                    var genomePos = ReferenceGenome.ConvertTorCRSPosition(position);
                    //insertion
                    if (refseq[i] == DnaAlphabet.Instance.Gap)
                    {
                        string insertedBases = "";
                        int j = i;
                        while (refseq[j] == DnaAlphabet.Instance.Gap)  {
                            if (qseq[j] == DnaAlphabet.Instance.Gap) { throw new HaploGrepException("Appears to be a shared gap in alignment, bad asssumption made!"); }
                            insertedBases += ((char)qseq[j]);
                            j++;
                        }
                        int insertionSize = j - i + 1;
                        //not sure whey the position is ever more than 1
                        string polyString=(genomePos-1).ToString()+".1"+insertedBases;
                        var p = new Polymorphism(polyString);
                        res.Add(p);
                        i = j - 1;
                    }
                    //deletion
                    else if (qseq[i] == DnaAlphabet.Instance.Gap)
                    {
                        //need to decide how big this guys is
                        int j = i;
                        while (qseq[j] == DnaAlphabet.Instance.Gap)
                        {
                            if (refseq[j] == DnaAlphabet.Instance.Gap) { throw new HaploGrepException("Appears to be a shared gap in alignment, bad asssumption made!"); }
                            string polyString = ReferenceGenome.ConvertTorCRSPosition(j).ToString() + "d";
                            var p = new Polymorphism(polyString);
                            res.Add(p);
                            j++; 
                        }
                        int deletionSize = j - i;
                        i=j-1;
                        position += deletionSize;
                    }
                    //SNP
                    else
                    {
                        var p = new Polymorphism(genomePos.ToString() + (char)qseq[i]);                        
                        res.Add(p);
                        position++;
                    }
                }
                else
                {
                    if (refseq[i] != DnaAlphabet.Instance.Gap) //should always be true
                    { position++; }
                }
            }
            return res;
        }
    
   
    }
}
