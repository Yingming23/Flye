import numpy 
import sys
with open("/raid/scratch/share/ont-haec/flye_example/features/deaf24fb-5ead-409a-87d4-468da8457089/chunk22_feats.npy", "rb") as f:
    a = numpy.load(f)
print(numpy.shape(a))
for i in range(len(a[0][0])):
    if a[1][0][i] != a[1][1][i]:
        print(i)

