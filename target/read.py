import sys
with open("target_726b2752-c218-4dbd-b3df-7d6733c58f4b.fastq", "r") as f:
    a = f.readlines()
seq = a[1]
chunk0 = seq[23137:36395]
with open("target_0914ddac-794f-41e9-bcf6-c3f9b4d0815b.fastq", "r") as f:
    b = f.readlines()
trg = b[1]
chunktrg = trg[94:13340]

with open("Analysis.txt", "w") as g:
    g.write(chunk0 + "\n")
    g.write(chunktrg + "\n")
    match = 0
    noMatch = ""
    non = False
    for j in range(13246):
        if chunk0[j] == chunktrg[j]:
            match += 1
        if chunk0[j] != chunktrg[j] and not non:
            noMatch += str(j) + " - "
            non = True
        elif chunk0[j] == chunktrg[j] and non:
            noMatch += str(j - 1) + ", "
            non = False
        elif j == 13245 and non:
            noMatch += str(j)
    g.write(str(match) + "\n")
    g.write(noMatch)
     
