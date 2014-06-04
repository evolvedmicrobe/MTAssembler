using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;
using Bio.IO;
using System.IO;
using Bio.IO.FastA;
using Bio.IO.FastQ;

namespace SimulateMitochondrialData
{
    class Program
    {
        static void Main(string[] args)
        {
            FastAParser fap = new FastAParser(@"D:\TestMixing\CRS.rn.fasta");
            FastQFormatter faq = new FastQFormatter(@"D:\TestMixing\RefSeq.fastq");
            var seq=fap.Parse().First() as Sequence;            
            seq=seq.GetSubSequence(1000, 300) as Sequence;
            int coverage = 50;
            int len = 75;
            byte[] QualScores = Enumerable.Range(0, len).Select(x => (byte)(33+35)).ToArray();
            for (int i = 0; i < (seq.Count - len); i++)
            {
                
                //string s=seq.GetSubSequence(i,len).ToString();
                string s=(seq as Sequence).ConvertToString(i,len );
                int hh = s.Length;
                byte[] sb= Encoding.UTF8.GetBytes(s);
                for (int j = 0; j < coverage; j++)
                { 
                    QualitativeSequence qs = new QualitativeSequence(Alphabets.DNA, FastQFormatType.Illumina_v1_8,sb, QualScores);
                    qs.ID="REF-"+i.ToString();
                    faq.Write(qs);
                }
            }
            faq.Close();
            fap.Close();
            byte[] newSeq=new byte[seq.Count];
            seq.CopyTo(newSeq,0,seq.Count);
            byte oldBase = seq[seq.Count / 2];
            newSeq[newSeq.Length/2]=65;
            //now mutate one base and go again
            Sequence nonRef=new Sequence(Alphabets.DNA,newSeq,true);
            faq = new FastQFormatter(@"D:\TestMixing\NonRefSeq.fastq");
            for (int i = 0; i < (seq.Count - len); i++)
            {   
                //string s=seq.GetSubSequence(i,len).ToString();
                string s=(nonRef as Sequence).ConvertToString(i,len);
                int hh = s.Length;
                byte[] sb= Encoding.UTF8.GetBytes(s);
                for (int j = 0; j < coverage; j++)
                { 
                    QualitativeSequence qs = new QualitativeSequence(Alphabets.DNA, FastQFormatType.Illumina_v1_8,sb, QualScores);
                    qs.ID="NONREF-"+i.ToString();
                    faq.Write(qs);
                }
            }
            faq.Close();
            
        }
    }
}
