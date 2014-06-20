from Bio import SeqIO;
filename=r"C:\Users\Delaney\SkyDrive\Software\MitochondrialPrograms\rCRS.fasta"
v=SeqIO.parse(open(filename),"fasta").next()
v=str(v.seq)
#should be 3106
v.index("N")