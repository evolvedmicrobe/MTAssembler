#r "C:\Program Files (x86)\.NET Bio\1.1\Tools\Bin\Bio.dll"
#r "C:\Program Files (x86)\.NET Bio\1.1\Tools\Bin\Bio.Pamsam.dll"

open Bio.IO.FastA
open Bio.Algorithms.Alignment.MultipleSequenceAlignment
open System
open System.IO
open Bio.SimilarityMatrices
open Microsoft.FSharp.Collections


let direc="D:\MTData\mtGenomes\DataFromPhyloTree\CombinedSequences\\"

let di=new DirectoryInfo(direc)
let files=di.GetFiles();
//Function to get files where the alphabet checks out ok
let getFasta (f:FileInfo) =
    try
        let parser=new FastAParser(f.FullName)
        parser.Alphabet<-Bio.Alphabets.DNA
        Some(parser.Parse() |> Seq.toArray)
    with
        | ex-> None
        
let q=getFasta(files.[0])

let FastaSeqs= files |> Array.Parallel.choose getFasta
                     |> Array.collect (fun x->x)
                     
let kmerLength=15
let distFunc=DistanceFunctionTypes.EuclideanDistance
let updater=UpdateDistanceMethodsTypes.Average
let profilerName=ProfileAlignerNames.SmithWatermanProfileAligner
let profileScore=ProfileScoreFunctionNames.InnerProduct
let similarityMatrix=new SimilarityMatrix(SimilarityMatrix.StandardSimilarityMatrix.AmbiguousDna)
let gapOpen= -5
let gapExtend= -2
let msa = new PAMSAMMultipleSequenceAligner(FastaSeqs,kmerLength,distFunc,
                                            updater,profilerName,profileScore,
                                            similarityMatrix,gapOpen,gapExtend,2,4)