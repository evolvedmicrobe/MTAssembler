using System;
using System.Collections;
using System.Collections.Generic;

namespace HaploGrep
{  
    public class Session
    {
        internal long timeLastAccess = 0L;
        internal SampleFile currentSampleFile;
        internal Dictionary<string, IList<ClusteredSearchResult>> classificationResults = new Dictionary<string, IList<ClusteredSearchResult>>();
        public Session()
        {
        }
        public virtual SampleFile CurrentSampleFile
        {
            get {
               return this.currentSampleFile;
            }
            set {
                this.currentSampleFile = value;
            }
        }
        public IList<ClusteredSearchResult> getClassificationResults(string sampleID)
        {
            return this.classificationResults[sampleID];
        }
        public void setClassificationResults(string sampleID, IList<ClusteredSearchResult> classificationResults)
        {
            this.classificationResults[sampleID] = classificationResults;

        }
        public virtual void resetClassificationResults()
        {
            this.classificationResults.Clear();
        }
        public void getUnusedPolys(string sampleID, string haplogroup)
        {
            throw new Exception("Need to possibly implement");
            //foreach (ClusteredSearchResult currentResult in (IList)this.classificationResults[sampleID])
            //{
                
            //    Element result = currentResult.getUnusedPolysXML(haplogroup);

            //    if (result != null)
            //    {
            //        return result;

            //    }
            //} return null;

        }


        //public virtual Element getNotInRangePolys(string sampleID, string haplogroup)
        //{
        //    foreach (ClusteredSearchResult currentResult in (IList)this.classificationResults[sampleID])
        //    {
        //        Element result = currentResult.getNotInRangePolysXML(haplogroup);

        //        if (result != null)
        //        {
        //            return result;

        //        }
        //    } return null;

        //}
    }

}