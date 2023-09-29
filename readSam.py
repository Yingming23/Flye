import pysam

pysam.sort("-o", "alignSort.sam", "alignment.sam")
samfile = pysam.AlignmentFile("alignSort.sam", "rb" )
samread = open("SamRead.txt", "w")
for pileupcolumn in samfile.pileup():
    samread.write("coverage at base %s = %s\n" % (pileupcolumn.pos, pileupcolumn.get_num_aligned()))
    for pileupread in pileupcolumn.pileups:
        # query position is None if is_del or is_refskip is set.
        samread.write('\tbase in read %s = %s\n' %
            (pileupread.alignment.query_name,
            pileupread.alignment.query_sequence[pileupread.query_position_or_next]))

samfile.close()
