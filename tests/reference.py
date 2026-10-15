"""Verbatim copies of the original implementations of functions that were optimized.

They are the specification the optimized versions are tested against (see test_*.py) and are not used by SplAdder itself.
"""

import scipy.sparse as spsp

from spladder.editgraph import *
from spladder.utils import *
from spladder.alt_splice.verify import *
from spladder.classes.counts import Counts


def insert_cassette_exons_reference(genes, bam_fnames, options):
    # written by Andre Kahles, Mpi Tuebingen, Germany, 2012

    inserted = 0

    strands = ['+', '-']

    ### form all possible combinations of contigs and strands --> regions
    (regions, options) = init_regions(bam_fnames, options.confidence, options, sparse_bam=options.sparse_bam)

    ### ignore contigs not present in bam files 
    # TODO TODO
    #keepidx = np.where(np.isin(np.array([CFG['chrm_lookup'][x.chr] for x in genes]), np.array([x.chr_num for x in regions])))[0]
    #genes = genes[keepidx]

    c = 0
    num_exons_added = 0
    num_exons = 0

    contigs = np.array([x.chr for x in genes], dtype='str')
    gene_strands = np.array([x.strand for x in genes])
    for contig in np.unique(contigs):
        bam_cache = dict()
        for si, s in enumerate(strands):
            cidx = np.where((contigs == contig) & (gene_strands == s))[0]

            for i in cidx:
                if options.verbose and (c+1) % 100 == 0:
                    print('\r %i(%i) genes done (found %i new cassette exons in %i tested intron pairs, %2.1f%%)' % (c+1, genes.shape[0], num_exons_added, num_exons, 100*num_exons_added/float(max(1, num_exons))))

                gg = genes[i]
                assert(gg.strand == s)
                assert(gg.chr == contig)

                if options.sparse_bam:
                    if isinstance(bam_fnames, str):
                        [tracks] = add_reads_from_sparse_bam(gg, bam_fnames, contig, options.confidence, types=['exon_track'], filter=options.read_filter, cache=bam_cache, unstranded=True)
                    else:
                        tracks = None
                        for fname in bam_fnames:
                            [tmp_] = add_reads_from_sparse_bam(gg, fname, contig, options.confidence, types=['exon_track'], filter=options.read_filter, cache=bam_cache, unstranded=True)
                            if tracks is None:
                                tracks = tmp_
                            else:
                                tracks += tmp_
                    tracks = np.asarray(tracks)
                else:
                    tracks = add_reads_from_bam(np.array([gg], dtype='object'), bam_fnames, ['exon_track'], options.read_filter, options.var_aware, options.primary_only, options.ignore_mismatches, mm_tag=options.mm_tag, cram_ref=options.ref_genome)

                ### add introns implied by splicegraph to the list
                all_introns = gg.introns[si][:, :2]
                for k in range(gg.splicegraph.edges.shape[0]):
                    for l in range(k+1, gg.splicegraph.edges.shape[0]):
                        if gg.splicegraph.edges[k, l] == 1:
                            all_introns = np.r_[all_introns, np.array([[gg.splicegraph.vertices[1, k], gg.splicegraph.vertices[0, l]]])] # introns are half open
                all_introns = unique_rows(all_introns)
           
                ### use only relevant introns (inside gene boundaries)
                if all_introns.shape[0] > 0:
                    keep_idx = np.where((all_introns[:, 1] > gg.start) & (all_introns[:, 0] < gg.stop))[0]
                    all_introns = all_introns[keep_idx, :]

                segment_starts = np.sort(np.unique(all_introns[:, 0]))
                segment_ends = np.sort(np.unique(all_introns[:, 1]))

                ### check for all intron-pairs, if exon could exist between them
                new_cassette = np.zeros((all_introns.shape[0], all_introns.shape[0]), dtype='int') 
                for k in range(all_introns.shape[0]):
                    for l in range(k + 1, all_introns.shape[0]):
                        if all_introns[k, 1] >= all_introns[l, 0]:
                            continue
                        ### only take intron pair, if outer ends are supported by current exons
                        if (not all_introns[k, 0] in gg.splicegraph.vertices[1, :]) or (not all_introns[l, 1] in gg.splicegraph.vertices[0, :]): 
                            continue
                        curr_exon = [all_introns[k, 1], all_introns[l, 0]]
                        ### do not allow curr_exon to overlap existing exon
                        if np.sum((gg.splicegraph.vertices[0, :] < curr_exon[1]) & (gg.splicegraph.vertices[1, :] > curr_exon[0])) > 0:
                            continue

                        num_exons += 1

                        if not ismember(curr_exon, gg.splicegraph.vertices.T, rows=True):
                            idx = np.arange(curr_exon[0], curr_exon[1]) - gg.start
                            exon_cov = np.sum(tracks[:, idx], axis=0)

                            pre_segment_end = np.where(segment_ends < curr_exon[0])[0]
                            if pre_segment_end.shape[0] > 0:
                                pre_segment_cov = np.sum(tracks[:, np.arange(segment_ends[pre_segment_end.max()], curr_exon[0]) - gg.start], axis=0)
                            else:
                                pre_segment_cov = np.sum(tracks[:, np.arange(curr_exon[0] - gg.start)], axis=0)
                            min_len_pre = min(pre_segment_cov.shape[0], exon_cov.shape[0])

                            aft_segment_start = np.where(segment_starts > curr_exon[1])[0]
                            if aft_segment_start.shape[0] > 0:
                                aft_segment_cov = np.sum(tracks[:, np.arange(curr_exon[1], segment_starts[aft_segment_start.min()]) - gg.start], axis=0)
                            else:
                                aft_segment_cov = np.sum(tracks[:, (curr_exon[1] - gg.start):], axis=0)
                            min_len_aft = min(aft_segment_cov.shape[0], exon_cov.shape[0])

                            if np.mean(exon_cov > (0.2 * np.mean(exon_cov))) > options.cassette_exon['min_cassette_region'] and \
                               np.median(exon_cov) > options.cassette_exon['min_cassette_cov']:
                                median_aft = max(np.median(aft_segment_cov[:min_len_aft]), 1)
                                median_pre = max(np.median(pre_segment_cov[-min_len_pre:]), 1)
                                if (np.median(exon_cov[-min_len_aft:]) / median_aft) - 1 >= options.cassette_exon['min_cassette_rel_diff'] and \
                                   (np.median(exon_cov[:min_len_pre]) / median_pre) - 1 >= options.cassette_exon['min_cassette_rel_diff']:
                                    new_cassette[k, l] = 1
                                    inserted += 1 
                any_added = False
                if any(new_cassette.ravel()):
                    curr_sg = gg.splicegraph.vertices
                    for k in range(new_cassette.shape[1]):
                        for l in range(k + 1, new_cassette.shape[1]):
                            if new_cassette[k, l] > 0:
                                exons_pre = np.where(curr_sg[1, :] == all_introns[k, 0])[0]
                                exons_aft = np.where(curr_sg[0, :] == all_introns[l, 1])[0]

                                gg.splicegraph.add_cassette_exon(np.array([all_introns[k, 1], all_introns[l, 0]]), exons_pre, exons_aft)
                                new_cassette[k, l] = 0
                                any_added = True
                                num_exons_added += 1
                    exon_order = np.argsort(gg.splicegraph.vertices[0, :])
                    gg.splicegraph.reorder(exon_order)
                    if not any_added:
                        break
                if any_added:
                    exon_order = np.argsort(gg.splicegraph.vertices[0, :])
                    gg.splicegraph.reorder(exon_order)

                ### clean up gene structure
                genes[i] = gg
                c += 1

    return (genes, inserted)


def label_alt_reference(self):

    # This script labels a gene as 'is_alt_spliced' that
    # has any transcript position that can be intronic or
    # exonic.
    #
    # This script labels each gene 'is_alt', that has a
    # simple alternative transcript start or end. Simple meaning
    # that the remaining transcript structure is the same and no 
    # other exons overlap to the alternative start or end.

    num_exons = self.splicegraph.vertices.shape[1]

    vertices = self.splicegraph.vertices
    edges = self.splicegraph.edges

    ### no edges in graph --> return
    if num_exons < 2:
        self.is_alt_spliced = 0
        self.is_alt = 0
        return

    init_alt_idx = []
    term_alt_idx = []
    for i in range(num_exons):
        ### handle start terminal exons -> no incoming edges
        if not np.any(edges[:i, i]):
            ### find other exons with same end
            idx = np.where(~np.isin(np.where(vertices[1, i] == vertices[1, :])[0], i))[0]
            if idx.shape[0] > 0:
                is_simple = True
                for j in idx:
                    ### not simple, if different introns follow --> break
                    if not isequal(np.where(edges[j, j+1:])[0] + j, np.where(edges[i, i+1:])[0] + i):
                        is_simple = False
                        break
                    ### not simple, if exons previous to j overlap the start of i --> break
                    idx_prev = np.where(edges[:j, j])[0]
                    for k in idx_prev:
                        if vertices[1, k] > vertices[0, i]:
                            is_simple = False
                            break
                if not is_simple:
                    continue
                ### if simple, report alternative initial exon
                init_alt_idx.append(i)
        ### handle end terminal exons -> no outgoing edges
        if not np.any(edges[i, i+1:]):
            ### find other exons with the same start
            idx = np.where(~np.isin(np.where(vertices[0, i] == vertices[0, :])[0], i))[0]
            if idx.shape[0] > 0: 
                is_simple = True
                for j in idx:
                    ### not simple, if different introns precede --> break
                    if not isequal(np.where(edges[:j, j])[0], np.where(edges[:i, i])[0]):
                        is_simple = False
                        break
                    ### not simple, if exons following to j start before i ends --> break
                    idx_next = np.where(edges[j, j+1:])[0] + j
                    for k in idx_next:
                        if vertices[0, k] < vertices[1, i]:
                            is_simple = False
                            break
                if not is_simple:
                    continue
                ### if simple, report alternative terminal exon
                term_alt_idx.append(i)

    ### further only consider exons that are neither init_alt nor term_alt
    take_idx = np.where(~np.isin(np.arange(num_exons), [init_alt_idx + term_alt_idx]))[0]
    if take_idx.shape[0] > 0:
        vertices = self.splicegraph.vertices[:, take_idx]
        edges = self.splicegraph.edges[take_idx, :][:, take_idx]

        start = vertices[0, :].min()
        stop = vertices[1, :].max()
        exon_loc = np.zeros((stop - start,))

        ### map all edges (introns) to genomic coordinates
        for i in range(edges.shape[0]):
            for j in range(i + 1, edges.shape[0]):
                if edges[i, j] == 1:
                    cur_edge = np.array([vertices[1, i], vertices[0, j]]) - start
                    exon_loc[cur_edge[0]:cur_edge[1]] += 1

        ### map all vertices (exons) to genomic coordinates 
        for i in range(vertices.shape[1]):
            cur_vertex = vertices[:, i] - start
            exon_loc[cur_vertex[0]:cur_vertex[1]] += 1
    else:
        exon_loc = np.zeros((1,))

    ### if at any position more than one exon or intron -> is_alt_spliced
    if max(exon_loc) > 1:
        self.is_alt_spliced = 1
        self.is_alt = 1
    else:
        self.is_alt_spliced = 0
        ### if not alt_spliced but term_alt or init_alt --> is_alt
        if len(init_alt_idx) == 0 and len(term_alt_idx) == 0:
            self.is_alt = 0 
        else:
            self.is_alt = 1


def verify_wrapper_reference(ev, genes, gidx_min, gene_ids_edges, gene_ids_segs, edge_idx, sample_idx, event_type, fn_count, options, idx):

    IN = h5py.File(fn_count, 'r')
    counts = []
    verified = []
    for i in range(ev.shape[0]):
        
        #sys.stdout.write('.')
        #if i > 0 and i % 50 == 0:
        #    sys.stdout.write('%i (%i)\n' % (i, ev.shape[0]))
        #sys.stdout.flush()
    
        g_idx = ev[i].gene_idx
        _g_idx = g_idx - gidx_min

        if genes[_g_idx].segmentgraph is None or genes[_g_idx].segmentgraph.segments.shape[1] == 0:
            genes[_g_idx].from_sparse()
            genes[_g_idx].segmentgraph = Segmentgraph(genes[_g_idx])
            genes[_g_idx].to_sparse()

        ### there are no edges present in the event
        if gene_ids_edges.shape[0] == 0:
            ver, info = verify_empty(event_type)
            counts.append(np.array([info]))
            verified.append([ver for _ in sample_idx])
            continue

        gr_idx_segs = np.where(gene_ids_segs == g_idx)[0]
        gr_idx_edges = np.where(gene_ids_edges == g_idx)[0]
        if gr_idx_edges.shape[0] == 0:
            ver, info = verify_empty(event_type)
            counts.append(np.array([info]))
            verified.append([ver for _ in sample_idx])
            continue

        if isinstance(sample_idx, int):
            sample_idx = [sample_idx]

        ### laod relevant count data from HDF5
        segments = np.atleast_2d(IN['segments'][gr_idx_segs, :])[:, sample_idx]
        seg_pos = np.atleast_2d(IN['seg_pos'][gr_idx_segs, :])[:, sample_idx]
        edges = np.atleast_2d(IN['edges'][gr_idx_edges, :])[:, sample_idx]
        curr_edge_idx = edge_idx[gr_idx_edges]

        verified.append([])
        for s_idx in range(len(sample_idx)):
            if event_type == 'exon_skip':
                ver, info = verify_exon_skip(ev[i], genes[_g_idx], segments[:, s_idx].T,  np.c_[curr_edge_idx, edges[:, s_idx]], options)
            elif event_type in ['alt_3prime', 'alt_5prime']:
                ver, info = verify_alt_prime(ev[i], genes[_g_idx], segments[:, s_idx].T,  np.c_[curr_edge_idx, edges[:, s_idx]], options)
            elif event_type == 'intron_retention':
                ver, info = verify_intron_retention(ev[i], genes[_g_idx], segments[:, s_idx].T,  np.c_[curr_edge_idx, edges[:, s_idx]], seg_pos[:, s_idx].T, options)
            elif event_type == 'mult_exon_skip':
                ver, info = verify_mult_exon_skip(ev[i], genes[_g_idx], segments[:, s_idx].T,  np.c_[curr_edge_idx, edges[:, s_idx]], options)
            elif event_type == 'mutex_exons':
                ver, info = verify_mutex_exons(ev[i], genes[_g_idx], segments[:, s_idx].T,  np.c_[curr_edge_idx, edges[:, s_idx]], options)

            verified[-1].append(ver)
            if s_idx == 0:
                counts.append(info[np.newaxis, :])
            else:
                counts[-1] = np.append(counts[-1], info[np.newaxis, :], axis=0) 
        verified[-1] = np.array(verified[-1], dtype='bool')

    IN.close()

    return (np.dstack(counts), np.dstack(verified), idx)


def segmentgraph_from_gene_reference(self, gene): 

    sg = gene.splicegraph.vertices
    breakpoints = np.unique(sg.ravel())
    self.segments = np.zeros((2, 0), dtype='int')
    for j in range(1, breakpoints.shape[0]):
        s = np.sum(sg[0, :] < breakpoints[j])
        e = np.sum(sg[1, :] < breakpoints[j])
        if s > e:
            self.segments = np.c_[self.segments, [breakpoints[j-1], breakpoints[j]]]

    ### match nodes to segments
    self.seg_match = np.zeros((0, sg.shape[1]), dtype='bool')
    for j in range(sg.shape[1]):
        tmp = ((sg[0, j] <= self.segments[0, :]) & (sg[1, j] >= self.segments[1, :]))
        if self.seg_match.shape[0] == 0:
            self.seg_match = tmp.copy().reshape((1, tmp.shape[0]))
        else:
            self.seg_match = np.r_[self.seg_match, tmp.reshape((1, tmp.shape[0]))]

    ### create edge graph between segments
    self.seg_edges = np.zeros((self.segments.shape[1], self.segments.shape[1]), dtype='bool')
    k, l = np.where(np.triu(gene.splicegraph.edges))

    for m in range(k.shape[0]):
        ### donor segment
        d = np.where(self.seg_match[k[m], :])[0][-1]
        ### acceptor segment
        a = np.where(self.seg_match[l[m], :])[0][0]
        self.seg_edges[d, a] = True


def count_gene_reference(gg, tracks, intron_list, intron_tol=0):
    """The original per gene counting of count.count_graph_coverage (verbatim loops), returns the Counts object."""

    counts = Counts(gg.segmentgraph.segments.shape[1])

    ### extract mean exon coverage for all segments
    for j in range(gg.segmentgraph.segments.shape[1]):
        idx = np.arange(gg.segmentgraph.segments[0, j], gg.segmentgraph.segments[1, j]) - gg.start
        counts.segments[j] = np.mean(np.sum(tracks[:, idx], axis=0))
        counts.seg_pos[j] = np.sum(np.sum(tracks[:, idx], axis=0) > 0)

    k, l = np.where(gg.segmentgraph.seg_edges == 1)

    ### there are no introns to count
    if intron_list.shape[0] == 0:
        for m in range(k.shape[0]):
            if counts.edges.shape[0] == 0:
                counts.edges = np.atleast_2d(np.array([np.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))
            else:
                counts.edges = np.r_[counts.edges, np.atleast_2d(np.array([np.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))]
        return counts

    ### extract intron counts 
    for m in range(k.shape[0]):
        idx = np.where((np.absolute(intron_list[:, 0] - gg.segmentgraph.segments[1, k[m]]) <= intron_tol) & (np.absolute(intron_list[:, 1] - gg.segmentgraph.segments[0, l[m]]) <= intron_tol))[0]
        if counts.edges.shape[0] == 0:
            if idx.shape[0] > 0:
                counts.edges = np.atleast_2d(np.array([np.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), np.sum(intron_list[idx, 2])]))
            else:
                counts.edges = np.atleast_2d(np.array([np.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))
        else:
            if idx.shape[0] > 0:
                counts.edges = np.r_[counts.edges, np.atleast_2d(np.array([np.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), np.sum(intron_list[idx, 2])]))]
            else:
                counts.edges = np.r_[counts.edges, np.atleast_2d(np.array([np.ravel_multi_index([k[m], l[m]], gg.segmentgraph.seg_edges.shape), 0]))]
    return counts


def to_sparse_reference(self):

    self.splicegraph_edges_shape = self.splicegraph.edges.shape
    self.splicegraph.edges = spsp.coo_matrix(self.splicegraph.edges) 
    self.splicegraph_edges_data = self.splicegraph.edges.data
    self.splicegraph_edges_row = self.splicegraph.edges.row
    self.splicegraph_edges_col = self.splicegraph.edges.col
    self.splicegraph.edges = None

    if hasattr(self, 'edge_count'):
        self.edge_count_shape = self.edge_count.shape
        self.edge_count = spsp.coo_matrix(self.edge_count)
        self.edge_count_data = self.edge_count.data
        self.edge_count_row = self.edge_count.row
        self.edge_count_col = self.edge_count.col
        self.edge_count = None 


def from_sparse_reference(self):

    if hasattr(self, 'splicegraph_edges_data'):
        self.splicegraph.edges = spsp.coo_matrix((self.splicegraph_edges_data, (self.splicegraph_edges_row, self.splicegraph_edges_col)), shape=self.splicegraph_edges_shape).toarray()
        del self.splicegraph_edges_data
        del self.splicegraph_edges_shape
        del self.splicegraph_edges_row
        del self.splicegraph_edges_col

    if hasattr(self, 'edge_count_data'):
        self.edge_count = spsp.coo_matrix((self.edge_count_data, (self.edge_count_row, self.edge_count_col)), shape=self.edge_count_shape).toarray()
        del self.edge_count_data
        del self.edge_count_shape
        del self.edge_count_row
        del self.edge_count_col


def insert_intron_retentions_reference(genes, bam_fnames, options):
    # written by Georg Zeller & Gunnar Raetsch, Mpi Tuebingen, Germany, 2009, Andre Kahles, MSKCC, 2013

    inserted = 0
    strands = ['+', '-']

    ### form all possible combinations of contigs and strands --> regions
    (regions, options) = init_regions(bam_fnames, options.confidence, options, sparse_bam=options.sparse_bam)

    ### ignore contigs not present in bam files 
    # TODO
    #keepidx = np.where(np.isin(np.array([options.chrm_lookup[x.chr] for x in genes[chunk_idx]]), np.array([x.chr_num for x in regions])))[0]
    #genes = genes[keepidx]

    c = 0 
    num_introns_added = 0
    num_introns = 0

    contigs = np.array([x.chr for x in genes], dtype='str')
    gene_strands = np.array([x.strand for x in genes])
    for contig in np.unique(contigs):
        bam_cache = dict()
        for si, s in enumerate(strands):
            cidx = np.where((contigs == contig) & (gene_strands == s))[0]

            for i in cidx:

                if options.verbose and (c+1) % 100 == 0:
                    print('\r %i(%i) genes done (found %i new retentions in %i tested eligible introns, %2.1f%%)' % (c+1, genes.shape[0], num_introns_added, num_introns, 100 * num_introns_added / float(max(1, num_introns))), file=sys.stdout)

                gg = genes[i]
                assert(gg.strand == s)
                assert(gg.chr == contig)

                if options.sparse_bam:
                    if isinstance(bam_fnames, str):
                        [tracks] = add_reads_from_sparse_bam(gg, bam_fnames, contig, options.confidence, types=['exon_track'], filter=options.read_filter, cache=bam_cache, unstranded=True)
                    else:
                        tracks = None
                        for fname in bam_fnames:
                            [tmp_] = add_reads_from_sparse_bam(gg, fname, contig, options.confidence, types=['exon_track'], filter=options.read_filter, cache=bam_cache, unstranded=True)
                            if tracks is None:
                                tracks = tmp_
                            else:
                                tracks += tmp_
                    tracks = np.asarray(tracks)
                else:
                    tracks = add_reads_from_bam(np.array([gg], dtype='object'), bam_fnames, ['exon_track'], options.read_filter, options.var_aware, options.primary_only, options.ignore_mismatches, mm_tag=options.mm_tag, cram_ref=options.ref_genome)

                exon_coverage = np.zeros((gg.splicegraph.vertices.shape[1],), dtype='float')
                for k in range(gg.splicegraph.vertices.shape[1]):
                    idx = np.arange(gg.splicegraph.vertices[0, k], gg.splicegraph.vertices[1, k]) - gg.start
                    exon_coverage[k] = np.median(np.sum(tracks[:, idx], axis=0).astype('float')) # median coverage for exon k

                ### check for all vertex-pairs, if respective intron can be retained
                new_retention = np.zeros(gg.splicegraph.edges.shape, dtype='int') 
                #for k in range(gg.splicegraph.edges.shape[0]):
                #    for l in range(k + 1, gg.splicegraph.edges.shape[0]):
                #        if gg.splicegraph.edges[k, l] == 1:
                for k,l in np.array(np.where(np.triu(gg.splicegraph.edges))).T:
                    ### ignore introns that include at least one complete exon
                    if np.sum((gg.splicegraph.vertices[0, :] > gg.splicegraph.vertices[1, k]) & (gg.splicegraph.vertices[1, :] < gg.splicegraph.vertices[0, l])) > 0:
                        continue

                    num_introns += 1
                    idx = np.arange(gg.splicegraph.vertices[1, k], gg.splicegraph.vertices[0, l]) - gg.start
                    icov = np.sum(tracks[:, idx], axis=0) 
                    if np.median(icov) > options.intron_retention['min_retention_cov'] and \
                        np.mean(icov > (0.5 * np.mean(icov))) > options.intron_retention['min_retention_region'] and  \
                        max(exon_coverage[k], exon_coverage[l]) / (1e-6 + min(exon_coverage[k], exon_coverage[l])) <= options.intron_retention['min_retention_max_exon_fold_diff'] and \
                        np.mean(icov) >= options.intron_retention['min_retention_rel_cov'] * (exon_coverage[k] + exon_coverage[l]) / 2.0 and \
                        np.mean(icov) <= options.intron_retention['max_retention_rel_cov'] * (exon_coverage[k] + exon_coverage[l]) / 2.0:

                        new_retention[k, l] = 1
                        inserted += 1
                    #	fprintf(log_fd, '%s\tintron_retention\t%c\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%2.1f\n', gg.chr, gg.strand, gg.splicegraph{1}(1,k), gg.splicegraph{1}(2,k), gg.splicegraph{1}(1,l), gg.splicegraph{1}(2,l), ...
                    #			floor(median(icov(1,:)+icov(2,:))), floor(gg.exon_coverage(k)), floor(gg.exon_coverage(l)), 100*mean(icov(1,:)+icov(2,:)>0)) ;
                any_added = False
                if False:
                    if np.sum(new_retention.ravel()) > 0:
                        new_retention = scipy.sparse.linalg.expm(new_retention)
                        new_retention[new_retention == 0] = 2
                        np.fill_diagonal(new_retention, 2)
                        while True:
                            any_added = False
                            k,l = np.unravel_index(new_retention.argmin(), new_retention.shape)
                            if new_retention[k, l] == 2:
                                break
                            if new_retention[k, l] > 0:
                                gg.splicegraph.add_intron_retention(k, l)
                                new_retention = np.c_[new_retention, np.ones((new_retention.shape[0], 1), dtype='int') * 2]
                                new_retention = np.r_[new_retention, np.ones((1, new_retention.shape[1]), dtype='int') * 2]
                                ### unset all inbetween retentions
                                for u in range(k, l + 1):
                                    for v in range(u + 1, l + 1):
                                        new_retention[u, v] = 2
                                #new_retention[k, l] = 0
                                any_added = True
                                num_introns_added += 1
                                #fprintf(log_fd, '%s\tintron_retention\t%i\t%i\t%i\t%i\t%i\t%2.1f\n', gg.chr, gg.splicegraph{1}(2,k), gg.splicegraph{1}(1,l), floor(median(icov(1,:)+icov(2,:))), gg.exon_coverage(k), gg.exon_coverage(l), 100*mean(icov(1,:)+icov(2,:)>0)) ;
                            exon_order = np.argsort(gg.splicegraph.vertices[0, :])
                            gg.splicegraph.reorder(exon_order)
                            new_retention = new_retention[exon_order, :][:, exon_order]
                            if not any_added:
                                break
                else:
                    if np.sum(new_retention.ravel()) > 0:
                        new_retention = scipy.sparse.linalg.expm(new_retention)
                        while True:
                            any_added = False
                            for k in range(new_retention.shape[1]):
                                for l in range(k + 1, new_retention.shape[1]):
                                    if new_retention[k, l] > 0:
                                        gg.splicegraph.add_intron_retention(k, l)
                                        new_retention = np.c_[new_retention, np.zeros((new_retention.shape[0], 1), dtype='int')]
                                        new_retention = np.r_[new_retention, np.zeros((1, new_retention.shape[1]), dtype='int')]
                                        new_retention[k, l] = 0
                                        any_added = True
                                        num_introns_added += 1
                                        #fprintf(log_fd, '%s\tintron_retention\t%i\t%i\t%i\t%i\t%i\t%2.1f\n', gg.chr, gg.splicegraph{1}(2,k), gg.splicegraph{1}(1,l), floor(median(icov(1,:)+icov(2,:))), gg.exon_coverage(k), gg.exon_coverage(l), 100*mean(icov(1,:)+icov(2,:)>0)) ;
                                        break
                                if any_added:
                                    break
                            exon_order = np.argsort(gg.splicegraph.vertices[0, :])
                            gg.splicegraph.reorder(exon_order)
                            new_retention = new_retention[exon_order, :][:, exon_order]
                            if not any_added:
                                break
                if any_added:
                    exon_order = np.argsort(gg.splicegraph.vertices[0, :])
                    gg.splicegraph.reorder(exon_order)
                genes[i] = gg
                c += 1
    return (genes, inserted)


def splicegraph_from_gene_reference(self, gene):
    
    for transcript_idx in range(len(gene.transcripts)):
        exon_start_end = gene.exons[transcript_idx]
        
        ### only one exon in the transcript
        if exon_start_end.shape[0] == 1:
            exon1_start = exon_start_end[0, 0]
            exon1_end = exon_start_end[0, 1]

            if self.vertices.shape[1] == 0:
                self.vertices = np.array([[exon1_start], [exon1_end]], dtype='int')
                self.edges = np.array([[0]], dtype='int')
            else:
                self.vertices = np.c_[self.vertices, [exon1_start, exon1_end]]
                self.new_edge()
        ### more than one exon in the transcript
        else:
            for exon_idx in range(exon_start_end.shape[0] - 1):
                exon1_start = exon_start_end[exon_idx , 0]
                exon1_end = exon_start_end[exon_idx, 1]
                exon2_start = exon_start_end[exon_idx + 1, 0]
                exon2_end = exon_start_end[exon_idx + 1, 1]
      
                if self.vertices.shape[1] == 0:
                    self.vertices = np.array([[exon1_start, exon2_start], [exon1_end, exon2_end]], dtype='int')
                    self.edges = np.array([[0, 1], [1, 0]], dtype='int')
                else:
                    exon1_idx = -1
                    exon2_idx = -1
                    ### check if current exon already occurred
                    for idx in range(self.vertices.shape[1]):
                        if ((self.vertices[0, idx] == exon1_start) and (self.vertices[1, idx] == exon1_end)):
                             exon1_idx = idx
                        if ((self.vertices[0, idx] == exon2_start) and (self.vertices[1, idx] == exon2_end)):
                             exon2_idx = idx

                    ### both exons already occured -> only add an edge
                    if (exon1_idx != -1) and (exon2_idx != -1):
                        self.edges[exon1_idx, exon2_idx] = 1
                        self.edges[exon2_idx, exon1_idx] = 1
                    else:
                        ### 2nd exon occured
                        if ((exon1_idx == -1) and (exon2_idx != -1)):
                            self.vertices = np.c_[self.vertices, [exon1_start, exon1_end]]
                            self.new_edge()
                            self.edges[exon2_idx, -1] = 1
                            self.edges[-1, exon2_idx] = 1
                        ### 1st exon occured
                        elif ((exon2_idx == -1) and (exon1_idx != -1)):
                            self.vertices = np.c_[self.vertices, [exon2_start, exon2_end]]
                            self.new_edge()
                            self.edges[exon1_idx, -1] = 1
                            self.edges[-1, exon1_idx] = 1
                        ### no exon occured
                        else:
                            assert((exon1_idx == -1) and (exon2_idx == -1))
                            self.vertices = np.c_[self.vertices, [exon1_start, exon1_end]]
                            self.vertices = np.c_[self.vertices, [exon2_start, exon2_end]]
                            self.new_edge()
                            self.new_edge()
                            self.edges[-2, -1] = 1
                            self.edges[-1, -2] = 1

    ### take care of the sorting by exon start
    s_idx = np.argsort(self.vertices[0, :])
    self.vertices = self.vertices[:, s_idx]
    self.edges = self.edges[s_idx, :][:, s_idx]
    self.terminals = np.zeros(self.vertices.shape, dtype='int')
    self.terminals[0, np.where(np.tril(self.edges).sum(axis=1) == 0)[0]] = 1
    self.terminals[1, np.where(np.triu(self.edges).sum(axis=1) == 0)[0]] = 1


def gene_add_exon_reference(self, exon, idx):
    if idx > (len(self.exons) - 1): 
        self.exons.append(np.zeros((0, 2), dtype='int'))
    self.exons[idx] = np.r_[self.exons[idx], np.expand_dims(exon, axis=0)]
