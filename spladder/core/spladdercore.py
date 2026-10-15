import sys
import os
import numpy as np
import pickle

from .. import evidence
from ..utils import dump_pickle, load_pickle
from .gen_graphs import gen_graphs 

def spladder_core(bam_fnames, out_fname, options):

    print('Augmenting splice graphs.', file=sys.stdout)
    print('=========================', file=sys.stdout)
    if not hasattr(options, 'genes'):
        genes = load_pickle(options.annotation, encoding='latin1')
    else:
        genes = options.genes

    ### mark which introns have been annotated
    for gene in genes:
        if not hasattr(gene, 'introns_anno'):
            gene.populate_annotated_introns()

    ### augment
    ### read evidence of the contigs is extracted in a single pass instead of once per gene and step
    with evidence.evidence_scope() as scope:
        evidence.prefetch_contigs(scope, bam_fnames, [gene.chr for gene in genes], options.read_filter, options, options.ignore_mismatches)
        genes = gen_graphs(genes, bam_fnames, options)

    print('Saving genes to %s' % (out_fname))
    dump_pickle(genes, out_fname)

