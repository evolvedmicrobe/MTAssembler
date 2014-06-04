import os
from Bio import SeqIO
refDirec=r"D:\MTData\MT_Data_Grabber\\"
refFile=r"human_g1k_v37.fasta"
queryFile=r"D:\MTData\MTNoN.fna"
os.chdir(refDirec)
infile=refFile
cmd="makeblastdb -in="+infile+" -input_type=fasta -dbtype=nucl"
print cmd
#os.system(cmd)
print "LLL"
cmd="blastn -query="+queryFile+" -outfmt=7 -out=BlastResult.txt -db="+infile
print(cmd)
#os.system(cmd)

outFile=open("BlastResult.txt")
v=[]
for line in outFile.readlines():
    isinstance(line,str)
    if line.startswith("#"):
        continue
    ls=line.split("\t")
    chrm=ls[1]
    if chrm=="MT":
        continue
    start=int(ls[8])
    end=int(ls[9])
    v.append(abs(end-start)+1)
    #print "chr"+chrm+":"+str(start)+"-"+str(end)
    print 'new MitoBlastMatchRegion("'+chrm+'",'+str(start)+","+str(end)+"),"
v.sort()
for i in v:
    print i