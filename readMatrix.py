import numpy 
import sys
with open("matrix1d92dd05-e7bb-407e-9428-1fd98dba6859.npy", "rb") as f:
    a = numpy.load(f)
print(numpy.shape(a))
with open("qry_id_contig_1.npy", "rb") as g:
    b = numpy.load(g)
chunk = 0
i = 0
while i != len(a[0][0]):
    total_len = 0
    cnt = 0
    while cnt != 4096:
        if total_len + i == len(a[0][0]):
            break
        if a[0][0][total_len + i] == ord("-"):
            cnt -= 1
        total_len += 1
        cnt += 1
    chunk_np = numpy.zeros((2, len(a[0]), total_len), dtype=numpy.uint8)
    for j in range(total_len):
        chunk_np[0][0][j] = a[0][0][i + j]
    aln_num = 1
    chunk_txt = open("chunk" + str(chunk) + "_qry_ids.txt", "w")
    for j in range(1, len(a[0])):
        if sum(a[0][j][i:total_len + i]) == 0:
            continue
        for k in range(total_len):
            chunk_np[0][aln_num][k] = a[0][j][k + i]
            chunk_np[1][aln_num][k] = a[1][j][k + i]
        for k in range(36):
            chunk_txt.write(b[j][k])
        chunk_txt.write("\n")
        aln_num += 1
    chunk_np = numpy.delete(chunk_np, numpy.s_[aln_num:], 1)
    with open("chunk" + str(chunk) + "_feats.npy", "wb") as chunk_arr:
        numpy.save(chunk_arr, chunk_np)
    chunk += 1
    i += total_len
    chunk_txt.close()

