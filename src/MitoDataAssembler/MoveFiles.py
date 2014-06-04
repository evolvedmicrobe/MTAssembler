import os
import shutil
h=os.getenv("HOME")
srcDir = os.path.join("bin","Release")
fileList=os.listdir(srcDir)
print fileList
fileList=[x for x in fileList if x.endswith(".dll") or x.endswith(".exe") or x.endswith(".pdb") or x.endswith(".mdb")]
for f in fileList:
	if f.count("vshost")>0:
		continue
	dest=os.path.join(h,"CSharpLib",f)
	src=os.path.join(srcDir,f)	
	shutil.copyfile(src,dest)
	print "Moved: "+f+" to cluster"
