from collections import namedtuple

# stores information for an allele, i.e., one variant
Allele = namedtuple("Allele",("loc","seq","freq","cov"))

class Haplotype(object):
    """Represents a unique haplotype found with some amount of read support.
    The read support may be poor, but will be reflected in the information
    contained by the object.  This will house the information to be used by
    the scoring algo.
    """
    def __init__(self):
        self.alleles = set()
        self.refid = 0
        self.refname = ''
        self.reads_expr = set()
        self.reads = set()
        self._prob = 0
        self._span = None

    def add_allele(self, allele):
        """Adds the given allele to the haplotype (if not already added)"""
        self.alleles.add(allele)

    def add_read(self, read_name, expr=False):
        """Add a read name to the haplotype. Recored any read that covers this
        haplotype, regardless of whether or not it's expressed within the read.
        If it is expressed, set expr=True.  Default is expr=False.
        """
        if expr:
            self.reads_expr.add(read_name)
        else:
            self.reads.add(read_name)
        
    def expressed(self, name):
        """Can only be used after the expressed reads have been recorded via th
        add_read method.  Returns true if the given read name is found in the
        expression set, otherwise false.
        """
        return name in self.reads_expr

    @property
    def all_reads(self):
        return self.reads | self.reads_expr

    @property
    def cov(self):
        return len(self.all_reads)

    @property
    def freq(self):
        return len(self.reads_expr)

    @property
    def prob(self):
        if not self._prob and self.cov:
            self._prob = self.freq / float(self.cov)

        return self._prob

    @property
    def span(self):
        """Returns a tuple (start, end) in genomic coordinates. Returns zero 
        if no alleles.
        """
        if len(self.alleles) == 0:
            return 0

        if not self._span:
            ordr = sorted(self.alleles,key=lambda a: a.loc)
            self._span = (ordr[0].loc, ordr[-1].loc)
        return self._span

    @property
    def hapstr(self):
        """String representation of the haplotype columns"""
        return "|".join(["%i-%s" % (a.loc, a.seq) 
                         for a in sorted(self.alleles,key=lambda a: a.loc)])

    def __repr__(self):
        astr = "%s freq: %i, cov: %i" % (self.hapstr, self.freq, self.cov)
        return astr

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.refid == other.refid and self.alleles == other.alleles
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(tuple([a for a in self.alleles]))
