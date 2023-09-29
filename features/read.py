import numpy
import sys
import string

an = open("Analysis.txt", "w")
flye = set()
ours = set()
M = set()
P = set()
N = set()
flyeLen = "("
flyeAln = "("
ourLen = "("
ourAln = "("
with open("/raid/scratch/share/ont-haec/yingming/chr19_m_read_ids.txt", "r") as aa:
    ab = aa.readlines()
    for i in ab:
        M.add(i[:36])
with open("/raid/scratch/share/ont-haec/yingming/chr19_p_read_ids.txt", "r") as ac:
    ad = ac.readlines()
    for i in ad:
        P.add(i[:36])
with open("/raid/scratch/liym/chr19_n_read_ids.txt", "r") as ae:
    af = ae.readlines()
    for i in af:
        N.add(i[:36])

for i in range(4):
    with open("chunk" + str(i) + "_feats.npy", "rb") as f:
        a = numpy.load(f)
    with open("chunk" + str(i) + "_ids.txt", "r") as x:
        c = x.readlines()
    with open("/raid/scratch/share/ont-haec/yingming/feats_mm2_acc/96d19b9b-a3ac-4202-95b0-1b4c635f5991/" + str(i) + ".ids.txt", "r") as y:
        d = y.readlines()
    flyeLen += str(numpy.shape(a)[1]) 
    flyeAln += str(len(c))
    if i != 3:
        flyeLen += ", "
        flyeAln += ", "
    else:
        flyeLen += ")"
        flyeAln += ")"
    with open("/raid/scratch/share/ont-haec/yingming/feats_mm2_acc/96d19b9b-a3ac-4202-95b0-1b4c635f5991/" + str(i) + ".features.npy", "rb") as g:
        b = numpy.load(g)
    ourLen += str(numpy.shape(b)[1])
    ourAln += str(len(d))
    if i != 3:
        ourLen += ", "
        ourAln += ", "
    else:
        ourLen += ")"
        ourAln += ")"
    for j in c:
        flye.add(j[:36])
    for j in d:
        ours.add(j[:36])
if "96d19b9b-a3ac-4202-95b0-1b4c635f5991" in M:
    an.write("M\n")
elif "96d19b9b-a3ac-4202-95b0-1b4c635f5991" in P:
    an.write("P\n")
else:
    an.write("N\n")
an.write("Flye:\n4 chunks of length " + flyeLen + "\n")
an.write(flyeAln + " alignments\n")
an.write("Ours:\n4 chunks of length " + ourLen + "\n")
an.write(ourAln + " alignments\n")
shared = set()
an.write("Alignments in common:\n")
for i in flye:
    if i in ours:
        shared.add(i)
for i in range(4):
    with open("chunk" + str(i) + "_feats.npy", "rb") as f:
        a = numpy.load(f)
    with open("/raid/scratch/share/ont-haec/yingming/feats_mm2_acc/96d19b9b-a3ac-4202-95b0-1b4c635f5991/" + str(i) + ".features.npy", "rb") as g:
        b = numpy.load(g)
    with open("chunk" + str(i) + "_ids.txt", "r") as x:
        c = x.readlines()
    with open("/raid/scratch/share/ont-haec/yingming/feats_mm2_acc/96d19b9b-a3ac-4202-95b0-1b4c635f5991/" + str(i) + ".ids.txt", "r") as y:
        d = y.readlines()
    an.write("Chunk " + str(i) + "\n")
    an.write("Flye:\n")
    for j in range(len(c)):
        if c[j][:36] in shared and j < 30:
            match = 0
            for k in range(len(a[0])):
                if chr(a[0][k][0]).lower() == chr(a[j + 1][k][0]).lower() and chr(a[0][k][0]).lower() != "-":
                    match += 1
            an.write("(Row " + str(j + 1) + ", ")
            if c[j][:36] in M:
                an.write("M, ")
            elif c[j][:36] in P:
                an.write("P, ")
            else:
                an.write("N, ")
                N.add(c[j][:36])
            an.write(c[j][:36].translate({ord(c): None for c in string.whitespace}) + ", Matched: " + str(match) + ")\n")
    an.write("Ours:\n")
    for j in range(len(d)):
        if d[j][:36] in shared and j < 30:
            match = 0
            for k in range(len(b[0])):
                if chr(b[0][k][0]).lower() == chr(b[j + 1][k][0]).lower() and chr(b[0][k][0]).lower() != "*":
                    match += 1
            an.write("(Row " + str(j + 1) + ", ")
            if d[j][:36] in M:
                an.write("M, ")
            elif d[j][:36] in P:
                an.write("P, ")
            else:
                an.write("N, ")
                N.add(d[j][:36])
            an.write(d[j].translate({ord(c): None for c in string.whitespace}) + ", Matched: " + str(match) + ")\n")
with open("/raid/scratch/liym/chr19_n_read_ids.txt", "w") as ae:
    for i in N:
        ae.write(i + "\n")

an.close()
