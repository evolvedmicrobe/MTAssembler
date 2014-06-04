This program is a set of code derived from Haplogrep that is designed to mimic the online tool.

The online tool was written in Java, and was somewhat complicated by the interface. In general, all the
complicated/historic methods are in the oldsearchmethods folder.  The NewSearchMethods folder contains
a more programmatic interface, a simpler rewrite, and loads the tree directly to avoid parsing the XML each 
time.

The Tests folder contains code to verify the mutations match the genbank accessions given by haplogrep (they don't).