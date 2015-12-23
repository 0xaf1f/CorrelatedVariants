def readMatch(variant, read):
    """True if the read contains the variant and the three flanking bases, 
    otherwise false.
    """
    vseq = variant.variantSeq
    vlen = len(vseq)
    algseq = read.read(orientation='genomic')
    refseq = read.reference(orientation='genomic')

    # locate the alignment offset, minding the genomic gaps
    offs = variant.start-(read.referenceStart+1)
    gaps = 0
    for idx, base in enumerate(refseq):
        if base == '-':
            gaps += 1
        if idx - gaps == offs:
            offs = idx
            break
            
    rloc = refseq[offs:offs+vlen]
    aloc = algseq[offs:offs+vlen]
    rflank = (refseq[offs-3:offs], refseq[offs+1:offs+4])
    aflank = (algseq[offs-3:offs], algseq[offs+1:offs+4])
    return vseq == aloc and aloc != rloc and rflank == aflank

def completeSpan(locus, reads):
    """Generator, returns reads that cover the given span completely."""
    for r in reads:
        if r[1] <= locus[1] and r[2] >= locus[2]:
            yield r
