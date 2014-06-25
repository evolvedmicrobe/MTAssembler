using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.Text;
using System.Linq;
using Bio.Util;

namespace Bio
{
    /// <summary>
    /// This class holds quality scores along with the sequence data.
    /// It is a container for the data in the SAMAlignedSequence, but much smaller to allow for faster parsing.
    /// much smaller.
    /// </summary>
    public class CompactSAMSequence : QualitativeSequence
    {
        #region SAMFields
        /// <summary>
        /// Name of Reference Sequence
        /// </summary>
        public string RName;
        /// <summary>
        /// Start position.
        /// </summary>
        public int Pos;
        /// <summary>
        /// SAM Flag vales
        /// </summary>
        public int SAMFlags;

        /// <summary>
        /// The SAM Cigar Field
        /// </summary>
        public string CIGAR;

        #endregion

        #region Constructors
        /// <summary>
        /// Initializes a new instance of the QualitativeSequence class with specified alphabet, quality score type,
        /// byte array representing symbols and encoded quality scores.
        /// Sequence and quality scores are validated with the specified alphabet and specified fastq format respectively.
        /// </summary>
        /// <param name="alphabet">Alphabet to which this instance should conform.</param>
        /// <param name="fastQFormatType">FastQ format type.</param>
        /// <param name="sequence">An array of bytes representing the symbols.</param>
        /// <param name="encodedQualityScores">An array of bytes representing the encoded quality scores.</param>
        public CompactSAMSequence(IAlphabet alphabet, FastQFormatType fastQFormatType, byte[] sequence, byte[] encodedQualityScores,bool validate)
            : base(alphabet, fastQFormatType, sequence, encodedQualityScores, validate)
        {
        }


        /// <summary>
        /// Initializes a new instance of the QualitativeSequence class with specified alphabet, quality score type,
        /// string representing symbols and encoded quality scores.
        /// Sequence and quality scores are validated with the specified alphabet and specified fastq format respectively.
        /// </summary>
        /// <param name="alphabet">Alphabet to which this instance should conform.</param>
        /// <param name="fastQFormatType">FastQ format type.</param>
        /// <param name="sequence">A string representing the symbols.</param>
        /// <param name="encodedQualityScores">A string representing the encoded quality scores.</param>
        public CompactSAMSequence(IAlphabet alphabet, FastQFormatType fastQFormatType, string sequence, string encodedQualityScores)
            : base(alphabet, fastQFormatType, sequence, encodedQualityScores, true)
        {
        }


        /// <summary>
        /// Initializes a new instance of the QualitativeSequence class with specified alphabet, quality score type,
        /// byte array representing symbols and signed byte array representing base quality scores 
        /// (Phred or Solexa base according to the FastQ format type).
        /// </summary>
        /// <param name="alphabet">Alphabet to which this instance should conform.</param>
        /// <param name="fastQFormatType">FastQ format type.</param>
        /// <param name="sequence">An array of bytes representing the symbols.</param>
        /// <param name="qualityScores">An array of signed bytes representing the base quality scores 
        /// (Phred or Solexa base according to the FastQ format type).</param>
        /// <param name="validate">If this flag is true then validation will be done to see whether the data is valid or not,
        /// else validation will be skipped.</param>
        public CompactSAMSequence(IAlphabet alphabet, FastQFormatType fastQFormatType, byte[] sequence, sbyte[] qualityScores, bool validate)
            : base(alphabet, fastQFormatType, sequence, qualityScores, validate) { }
        /// <summary>
        /// Initializes a new instance of the QualitativeSequence class with specified alphabet, quality score type,
        /// byte array representing symbols and integer array representing base quality scores 
        /// (Phred or Solexa base according to the FastQ format type).
        /// </summary>
        /// <param name="alphabet">Alphabet to which this instance should conform.</param>
        /// <param name="fastQFormatType">FastQ format type.</param>
        /// <param name="sequence">An array of bytes representing the symbols.</param>
        /// <param name="qualityScores">An array of integers representing the base quality scores 
        /// (Phred or Solexa base according to the FastQ format type).</param>
        /// <param name="validate">If this flag is true then validation will be done to see whether the data is valid or not,
        /// else validation will be skipped.</param>
        public CompactSAMSequence(IAlphabet alphabet, FastQFormatType fastQFormatType, byte[] sequence, int[] qualityScores, bool validate)
            : base(alphabet, fastQFormatType, sequence, qualityScores, validate) { }

        #endregion


    }
}
