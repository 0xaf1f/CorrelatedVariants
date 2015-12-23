# Scoring algorithm for the correlation.
import logging
import math
import itertools
from itertools import ifilter
import numpy
from scipy.stats import binom

class ScoreBase(object):
    """Abstract base-class defining the basic API contract for scoring a set of 
    potentially correlated variants. Mutual information can be computed many 
    different ways, so we generalize the API here.  Subclasses will explore 
    the various way to compute and represent the value.
    """
    def __init__(self):
        self.hapmap = dict()
        self.results = []
        self.score = False

    def correlate(self, score=False):
        """Builds the correlation data for all haplotypes. Don't score the
        correlation unless asked to (score=True). 
        """
        self.score = score
        for hap in ifilter(lambda h: len(h.alleles) > 1 and h.prob > 0.01 
                           ,self.hapmap.values()):
            line = [hap.refname, hap.hapstr, str(hap.freq), str(hap.cov), 
                    '%.2f' % (hap.freq/float(hap.cov)*100)]
            if self.score:
                line.extend(self.calculate(hap))

            self.results.append(line)

    def write(self, out):
        """Writes out the results in CSV format to the given output handle."""
        headers = ['ref','haplotype','frequency','coverage','percent']
        if self.score:
            headers.append('mutinf')
            headers.append('pval')

        out.write(",".join(headers))
        out.write("\n")

        # output sorted descending by haplotype percent
        for r in sorted(self.results, key=lambda r: float(r[4]), reverse=True):
            out.write(",".join(str(v) for v in r) + "\n")

    def calculate(self, haplotype):
        """The meat of the scoring algorithm. Left to implementors who must
        document their respective methods.  This function must return a tuple
        of mutual information score for the given haplotype, and a p-value
        indicating how likely this relationship is arises by chance.
        """
        pass

class TotalCorrelation(ScoreBase):
    """Multivariate generalization of the mutual information calculation.
    aka, multivariant constraint or multiinformation. Scales to any number
    of variants.  Again, the sets we're working with are {base, not base}

        C(X_1,X_2,\ldots,X_n) = \sum_{x_1\in\mathcal{X}_1}
            \sum_{x_2\in\mathcal{X}_2} \ldots \sum_{x_n\in\mathcal{X}_n}
            p(x_1,x_2,\ldots,x_n)\log\frac{p(x_1,x_2,\ldots,x_n)}
            {p(x_1)p(x_2)\cdots p(x_n)}

    """
        
    def calculate(self, haplotype):
        """ Returns: (mutinf, pval)
        """
        ndims = len(haplotype.alleles)
        hapstr = haplotype.hapstr

        # N-dimensional contingency table, one dimension for each allele
        # holds the observed frequencies.
        obs = numpy.zeros([2] * ndims)

        # Holds the expected frequencies, or product of individual frequencies.
        exp = numpy.zeros([2] * ndims)

        # list of possible alleles
        hsa = numpy.array(haplotype.hapstr.split("|"))

        # Iterate over the truth table and assign the probabilities.
        # Each row corresponds to a cell in the contigency table as
        # well as a possible haplotype combination.  Probabilities
        # must sum to 1.
        rem = 1
        perobs = self._permute(ndims)
        for row in perobs:
            # lookup string for haplotype combos
            lus = "|".join(hsa[numpy.array(row,dtype=bool)])

            # special case, no alleles present. compute later.
            if not lus:
                continue

            # check the hapmap for the observed probability, any missing
            # observations will drain into the remainder.
            if lus in self.hapmap:
                sub_hap = self.hapmap[lus]
                freq = len(haplotype.reads_expr & sub_hap.reads_expr)
                cov = float(len(haplotype.all_reads & sub_hap.all_reads))
                prob = freq / cov
                obs[row] = prob
                rem -= prob

        # remainder holds all negative observations
        # Bug 21042: careful, this went negative
        obs[tuple([0]*ndims)] = rem 

        # compute margins
        margins = self._margins(obs)

        # compute total correlation
        jiprs = []
        for row in perobs:
            # joint probablity
            jp = obs[row]
            # individual probablity product (numpy indexing is cool ...)
            ipp = numpy.prod(numpy.diagonal(margins[:,row]))
            exp[row] = ipp
            # i.e., observed
            if jp > 0:
                jiprs.append(jp*math.log(jp/ipp,2))

        cov = haplotype.cov
        row = tuple([1]*ndims)
        nsucc = obs[row]*cov
        hprob = numpy.prod(numpy.diagonal(margins[:,row]))
        
        # return (score, probability)
        return (sum(jiprs), binom.sf(nsucc, cov, hprob))

    def _margins(self, a):
        """0.9 does not have scipy.stats.contingency. Borrows from the scipy 
        0.11 stats.contingency package
        """
        margsums = numpy.zeros((a.ndim,2))
        _range = range(a.ndim)
        for k in xrange(a.ndim):
            marg = numpy.apply_over_axes(numpy.sum, a, [j for j in _range if j != k])
            margsums[k] = numpy.ndarray.flatten(marg)
        return margsums

    def _permute(self, ndims):
        """Generates a set of tuples containing all possible permutations of 
        the presence/absence of a given base in the haplotype. 
        """
        cols = [0]*ndims
        tt = set()
        for i in xrange(ndims):
            for row in itertools.permutations(cols):
                tt.add(row)
            cols[i] = 1

        tt.add(tuple([1]*ndims))
        return tt
