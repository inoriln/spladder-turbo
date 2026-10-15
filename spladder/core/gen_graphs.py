import numpy as np
import sys
import copy

if __package__ is None:
    __package__ = 'modules.core'

from .. import parallel
from ..init import init_regions
from ..reads import *
from ..helpers import *
from ..editgraph import *
from ..merge import *

def _enter_filter_swap(options, settings):
    """Some steps use read filters of their own: make them the current ones; returns what _leave_filter_swap needs."""

    CFG_ = dict()
    if 'read_filter' in settings:
        CFG_['read_filter'] = options.read_filter.copy()
        options.read_filter = settings['read_filter']
    return CFG_


def _leave_filter_swap(options, CFG_):
    for key in CFG_:
        setattr(options, key, CFG_[key].copy())


def _augment_genes(genes, introns, bam_fnames, options, inserted, fd_log):
    """Inserts cassette exons, intron retentions and new intron edges into the graphs of the genes and labels them again.

    The graph of a gene is edited with the read evidence of the gene only, so this can be done for any set of
    genes, e.g., for the genes of one contig (see _augment_genes_in_shards).  `introns` are the introns of the
    genes (None if none of the steps needs them), `inserted` counts what has been inserted."""

    ### insert casssette exons
    if options.insert_es:
        print('Inserting cassette exons ...', file=fd_log)
        CFG_ = _enter_filter_swap(options, getattr(options, 'cassette_exon', dict()))
        genes, inserted_ = insert_cassette_exons(genes, bam_fnames, options)
        inserted['cassette_exon'] = inserted_
        _leave_filter_swap(options, CFG_)
        print('\n... inserted %i casette exons ....\n... done.\n' % inserted['cassette_exon'], file=fd_log)

    if options.insert_ir:
        print('Inserting intron retentions ...', file=fd_log)
        CFG_ = _enter_filter_swap(options, options.intron_retention)
        genes, inserted_ = insert_intron_retentions(genes, bam_fnames, options)
        inserted['intron_retention'] = inserted_
        _leave_filter_swap(options, CFG_)
        print('\n... inserted %i new intron retentions ...\n...done.\n' % inserted['intron_retention'], file=fd_log)

    if options.remove_se:
        print('Removing short exons ...', file=fd_log)
        genes = remove_short_exons(genes, options)
        for i in range(genes.shape[0]):
            if np.any(genes[i].splicegraph.vertices[1, :] - genes[i].splicegraph.vertices[0, :] < options.remove_exons['min_exon_len_remove']):
                print('WARNING: could not remove all short exons', file=sys.stderr)
        print('... done.\n', file=fd_log)


    # test all exons if the reading frame is larger if exon is skipped
    ##############################################################################%%
    #print('find exons to skip to elongate reading frame', file=fd_log)
    #genes = insert_cds_exon_skips(genes, genome_info)
    #genes = splice_graph(genes)


    # sanity checking
    for g in genes:
        assert(all(g.splicegraph.vertices[0, :] <= g.splicegraph.vertices[1, :]))

    ### insert new intron edges
    if options.insert_ni:
        # re-set list of introns supported by RNA-seq data to 
        # the genes structure
        ##############################################################################%%
        for i in range(genes.shape[0]):
            genes[i].introns = introns[i, :]

        print('Inserting new intron edges ...', file=fd_log)
        chrms = np.array([x.chr for x in genes], dtype='str')
        for chr_idx in np.unique(chrms):
            genes_before = genes[np.where(chrms == chr_idx)[0]]
            tmp_genes = copy.deepcopy(genes_before)
            #
            ##############################################################################%%
            if not hasattr(options, 'insert_intron_iterations'):
                options.insert_intron_iterations = 5
            for iter in range(1, options.insert_intron_iterations + 1):
                print('... chr %s - iteration %i/%i\n' % (chr_idx, iter, options.insert_intron_iterations), file=fd_log)

                genes_mod, inserted_ = insert_intron_edges(tmp_genes, bam_fnames, options)

                inserted['intron_in_exon'] += inserted_['intron_in_exon']
                inserted['alt_53_prime'] += inserted_['alt_53_prime']
                inserted['exon_skip'] += inserted_['exon_skip']
                inserted['gene_merge'] += inserted_['gene_merge']
                inserted['new_terminal_exon'] += inserted_['new_terminal_exon']

                # in case any exon was inserted that already existed, we merge them into one exon 
                print('... removing duplicate exons ...', file=fd_log)
                genes_mod = merge_duplicate_exons(genes_mod, options)

                # inserted
                if isequal(genes_mod, genes_before):
                    break
                else:
                    _updated = np.where([genes_mod[i] != genes_before[i] for i in range(len(genes_mod))])[0].shape[0]
                    print('... updated %i genes' % _updated)

                genes_before = copy.deepcopy(genes_mod)
                tmp_genes = genes_mod
            chrms = np.array([x.chr for x in genes], dtype='str')
            genes[np.where(chrms == chr_idx)[0]] = copy.deepcopy(genes_mod)
        print('... done.\n', file=fd_log)

    print('Re-labeling new alternative genes ...', file=fd_log)
    for ix in range(genes.shape[0]):
        genes[ix].start = genes[ix].splicegraph.vertices.min()
        genes[ix].stop = genes[ix].splicegraph.vertices.max()
        genes[ix].label_alt()
    print('... done.\n', file=fd_log)

    return genes


def _augment_shard(task):
    """_augment_genes for the genes of one shard (the work of a worker process of _augment_genes_in_shards)."""

    genes, introns, bam_fnames, options, positions, keys = parallel.shared()
    pos = positions[task]
    shard_genes = genes[pos]
    shard_introns = None if introns is None else introns[pos, :]
    shard_inserted = dict((key, 0) for key in keys)
    fd_log = sys.stdout if options.logfile == '-' else open('%s.shard%i' % (options.logfile, task), 'w')
    try:
        shard_genes = _augment_genes(shard_genes, shard_introns, bam_fnames, options, shard_inserted, fd_log)
    finally:
        if fd_log is not sys.stdout:
            fd_log.close()
    return (pos, shard_genes, shard_inserted)


def _augment_genes_in_shards(genes, introns, bam_fnames, options, inserted, fd_log):
    """_augment_genes with the genes of different contigs being edited in different processes.

    The edits of a gene depend on the genes of its contig at most (the genes next to it), and the result does not
    depend on the order in which the contigs are processed.  What the steps leave behind in the settings is done here."""

    chrms = np.array([x.chr for x in genes], dtype='str')
    contigs = [(chrm, np.where(chrms == chrm)[0]) for chrm in np.unique(chrms)]
    workers = min(options.parallel, len(contigs))
    ### the contigs are put together to about as many shards as are good for balancing the load, the biggest first
    positions = parallel.balanced_shards(contigs, workers)

    print('Editing the splice graphs of %i contigs in %i processes ...' % (len(contigs), workers), file=fd_log)
    for (pos, shard_genes, shard_inserted) in parallel.run_shards(_augment_shard, (genes, introns, bam_fnames, options, positions, list(inserted.keys())), len(positions), workers):
        genes[pos] = shard_genes
        for key in inserted:
            inserted[key] += shard_inserted[key]

    ### the steps change the settings they use (see _augment_genes), and let the contigs of the alignment file be known
    if options.insert_es or options.insert_ir:
        init_regions(bam_fnames, options.confidence, options, sparse_bam=options.sparse_bam)
    if options.insert_es:
        _leave_filter_swap(options, _enter_filter_swap(options, getattr(options, 'cassette_exon', dict())))
        print('... inserted %i casette exons ....\n... done.\n' % inserted['cassette_exon'], file=fd_log)
    if options.insert_ir:
        _leave_filter_swap(options, _enter_filter_swap(options, options.intron_retention))
        print('... inserted %i new intron retentions ...\n...done.\n' % inserted['intron_retention'], file=fd_log)
    if options.insert_ni:
        if not hasattr(options, 'debug'):
            options.debug = False
        if not hasattr(options, 'insert_intron_iterations'):
            options.insert_intron_iterations = 5
    print('... done.\n', file=fd_log)
    return genes


def _can_edit_in_shards(genes, options):
    """The graphs of the contigs are edited in different processes, if there are processes to do so and nothing couples the genes of different contigs."""

    return (options.parallel > 1 and parallel.can_fork() and not options.remove_se and not options.sparse_bam and
            len(set(x.chr for x in genes)) > 1)


def gen_graphs(genes, bam_fnames, options=None):

    if options.logfile == '-':
        fd_log = sys.stdout
    else:
        fd_log = open(options.logfile, 'w')

    ### init the stats for inserted elements
    inserted = dict()
    inserted['cassette_exon'] = 0
    inserted['intron_retention'] = 0
    inserted['intron_in_exon'] = 0
    inserted['alt_53_prime'] = 0
    inserted['exon_skip'] = 0
    inserted['gene_merge'] = 0
    inserted['new_terminal_exon'] = 0

    # build splice graph for all genes 
    ##############################################################################%%
    print('Generating splice graph ...', file=fd_log)
    ### merge exons if possible / reduce graph
    ### originially implemented for ESTs, reduces complexity, 
    ### but removes alternative transcript starts and ends !
    if options.infer_sg:
        genes = infer_splice_graph(genes)

    ### sort exons by start position in ascending order
    for ix in range(genes.shape[0]):
        genes[ix].splicegraph.sort()

    ### label alternative and constitutive genes
    for ix in range(genes.shape[0]):
        genes[ix].label_alt()
    if options.verbose:
        print('\nTotal genes:\t\t\t\t\t\t\t%d' % genes.shape[0], file=fd_log)
        print('Total genes with alternative isoforms:\t\t\t\t%d' % np.sum([x.is_alt for x in genes]), file=fd_log)
        print('Total genes alternatively spliced:\t\t\t\t%d' % np.sum([x.is_alt_spliced for x in genes]), file=fd_log)
        print('Total constitutively spliced:\t\t\t\t\t%d' % (genes.shape[0] - np.sum([x.is_alt_spliced for x in genes])), file=fd_log)

    ### update terminals, start terminals in row 1, end terminals in row 2
    for i in range(genes.shape[0]):
        genes[i].splicegraph.update_terminals()

    ### reset gene start and gene stop according to exons and splice graph
    for i in range(genes.shape[0]):
        genes[i].start = min([x.min() for x in genes[i].exons])
        genes[i].stop = max([x.max() for x in genes[i].exons])
    print('...done.\n', file=fd_log)

    ### sort genes by positions
    s_idx = np.argsort([x.stop for x in genes])
    genes = genes[s_idx]
    s_idx = np.argsort([x.start for x in genes], kind='mergesort')
    genes = genes[s_idx]
    s_idx = np.argsort([x.chr for x in genes], kind='mergesort')
    genes = genes[s_idx]

    # append list of introns supported by RNA-seq data to 
    # the genes structure (only in case we want to augment the graph)
    ##############################################################################
    if (options.insert_es or options.insert_ir or options.insert_ni): 
        print('Loading introns from file ...', file=fd_log)
        introns = get_intron_list(genes, bam_fnames, options)
        print('...done.\n', file=fd_log)

        ### clean intron list
        ### remove all introns that do not conform to the splice site consensus
        if len(options.filter_consensus) > 0:
            print('Filtering introns for adherence to %s splice site consensus ...' % options.filter_consensus, file=fd_log)
            introns = filter_introns_consensus(introns, genes, options)
            print('...done.\n', file=fd_log)

        ### remove all introns that overlap more than one gene on the same strand
        print('Filtering introns for ambiguity ...', file=fd_log)
        introns = filter_introns(introns, genes, options)
        print('...done.\n', file=fd_log)

        ### check feasibility
        ### TODO when working exclusively with sparse bam, we need to skip this ...
        print('Testing for infeasible genes ...', file=fd_log)
        introns = make_introns_feasible(introns, genes, bam_fnames, options)
        print('...done.\n', file=fd_log)

        for i in range(genes.shape[0]):
            genes[i].introns = introns[i, :]

    introns = introns if (options.insert_es or options.insert_ir or options.insert_ni) else None
    if _can_edit_in_shards(genes, options):
        genes = _augment_genes_in_shards(genes, introns, bam_fnames, options, inserted, fd_log)
    else:
        genes = _augment_genes(genes, introns, bam_fnames, options, inserted, fd_log)

    ### print summary to log file
    print('Inserted:', file=fd_log)
    for fn in inserted:
        print('\t%s:\t%i' % (fn, inserted[fn]), file=fd_log)

    if fd_log != sys.stdout:
        fd_log.close()

    return (genes, inserted)

