import numpy as np
import scipy.sparse as spsp
import warnings

if __package__ is None:
    __package__ = 'modules.classes'

from ..utils import *
from .segmentgraph import Segmentgraph
from .splicegraph import Splicegraph

def _matrix_to_coo(matrix):
    """The non-zero entries of a matrix as (values, rows, columns) in row-major order, which is what scipy.sparse.coo_matrix(matrix) stores.

    The adjacency matrices of the splice graphs are plain two-dimensional arrays; this is much faster than going through scipy."""

    if type(matrix) is np.ndarray and matrix.ndim == 2 and matrix.dtype.kind in 'biuf' and max(matrix.shape) < 2 ** 31:
        row, col = matrix.nonzero()
        return (matrix[row, col], row.astype(np.int32), col.astype(np.int32))
    coo = spsp.coo_matrix(matrix)
    return (coo.data, coo.row, coo.col)


def _matrix_from_coo(data, row, col, shape):
    """The dense matrix of the entries stored by _matrix_to_coo."""

    if type(data) is np.ndarray and data.ndim == 1 and data.dtype.kind in 'biuf' and len(shape) == 2:
        matrix = np.zeros(shape, dtype=data.dtype)
        matrix[row, col] = data
        return matrix
    return spsp.coo_matrix((data, (row, col)), shape=shape).toarray()


### graphs with more exons than this are labelled by the general implementation (the compiled one compares all pairs of exons)
_MAX_COMPILED_EXONS = 2000


@jit_cached
def _same_successors(edges, a, b, n):
    """True, if the exons a and b have the same following exons (the columns of the non-zero entries right of the diagonal)."""

    ca = a + 1
    cb = b + 1
    while True:
        while ca < n and edges[a, ca] == 0:
            ca += 1
        while cb < n and edges[b, cb] == 0:
            cb += 1
        if ca >= n or cb >= n:
            return ca >= n and cb >= n
        if ca != cb:
            return False
        ca += 1
        cb += 1


@jit_cached
def _same_predecessors(edges, a, b):
    """True, if the exons a and b have the same preceding exons (the rows of the non-zero entries above the diagonal)."""

    ra = 0
    rb = 0
    while True:
        while ra < a and edges[ra, a] == 0:
            ra += 1
        while rb < b and edges[rb, b] == 0:
            rb += 1
        if ra >= a or rb >= b:
            return ra >= a and rb >= b
        if ra != rb:
            return False
        ra += 1
        rb += 1


@jit_cached
def _label_alt_kernel(vertices, edges):
    """(status, is_alt_spliced, is_alt) of a splice graph; the same logic as Gene._label_alt_numpy, step by step.

    status 1: the graph has no position covered by an exon (the general implementation fails for it, which the caller repeats),
    status 2: an exon ends before it starts."""

    n = vertices.shape[1]

    ### exons that end before they start: left to the general implementation
    for m in range(n):
        if vertices[1, m] < vertices[0, m]:
            return (2, 0, 0)

    ### only the upper triangle of the adjacency matrix is looked at
    has_in = np.zeros(n, dtype=np.bool_)
    has_out = np.zeros(n, dtype=np.bool_)
    shares_end = np.zeros(n, dtype=np.bool_)
    shares_start = np.zeros(n, dtype=np.bool_)
    for i in range(n):
        for j in range(i + 1, n):
            if edges[i, j] != 0:
                has_out[i] = True
                has_in[j] = True
            if vertices[1, i] == vertices[1, j]:
                shares_end[i] = True
                shares_end[j] = True
            if vertices[0, i] == vertices[0, j]:
                shares_start[i] = True
                shares_start[j] = True

    excluded = np.zeros(n, dtype=np.bool_)
    n_alt = 0
    for i in range(n):
        if not ((not has_in[i] and shares_end[i]) or (not has_out[i] and shares_start[i])):
            continue

        ### handle start terminal exons -> no incoming edges
        incoming = False
        for k in range(i):
            if edges[k, i] != 0:
                incoming = True
                break
        if not incoming:
            ### the other exons with the same end; they are numbered by their position in the list of the exons with that end (sic)
            has_other = False
            is_simple = True
            pos = 0
            for m in range(n):
                if vertices[1, m] == vertices[1, i]:
                    if m != i:
                        has_other = True
                        j = pos
                        ### not simple, if different introns follow
                        if not _same_successors(edges, j, i, n):
                            is_simple = False
                            break
                        ### not simple, if exons previous to j overlap the start of i
                        for k in range(j):
                            if edges[k, j] != 0 and vertices[1, k] > vertices[0, i]:
                                is_simple = False
                                break
                        if not is_simple:
                            break
                    pos += 1
            if has_other:
                if not is_simple:
                    continue
                excluded[i] = True
                n_alt += 1

        ### handle end terminal exons -> no outgoing edges
        outgoing = False
        for k in range(i + 1, n):
            if edges[i, k] != 0:
                outgoing = True
                break
        if not outgoing:
            has_other = False
            is_simple = True
            pos = 0
            for m in range(n):
                if vertices[0, m] == vertices[0, i]:
                    if m != i:
                        has_other = True
                        j = pos
                        ### not simple, if different introns precede
                        if not _same_predecessors(edges, j, i):
                            is_simple = False
                            break
                        ### not simple, if exons following to j start before i ends; the following exons are numbered one too low (sic)
                        for k in range(j + 1, n):
                            if edges[j, k] != 0 and vertices[0, k - 1] < vertices[1, i]:
                                is_simple = False
                                break
                        if not is_simple:
                            break
                    pos += 1
            if has_other:
                if not is_simple:
                    continue
                excluded[i] = True
                n_alt += 1

    ### further only consider exons that are neither init_alt nor term_alt
    n_take = 0
    lo = 0
    hi = 0
    for m in range(n):
        if not excluded[m]:
            if n_take == 0:
                lo = vertices[0, m]
                hi = vertices[1, m]
            else:
                lo = min(lo, vertices[0, m])
                hi = max(hi, vertices[1, m])
            n_take += 1

    overlapping = False
    if n_take > 0:
        if hi <= lo:
            return (1, 0, 0)
        ### Is any position covered by more than one intron or exon?  Introns are the stretches between the end of an
        ### exon and the start of a following exon it is connected to: look for two (half open) intervals that overlap.
        n_edges = 0
        for a in range(n):
            if excluded[a]:
                continue
            for b in range(a + 1, n):
                if not excluded[b] and edges[a, b] == 1:
                    n_edges += 1
        starts = np.empty(n_edges + n_take, dtype=np.int64)
        stops = np.empty(n_edges + n_take, dtype=np.int64)
        cnt = 0
        for a in range(n):
            if excluded[a]:
                continue
            for b in range(a + 1, n):
                if not excluded[b] and edges[a, b] == 1:
                    if vertices[0, b] > vertices[1, a]:
                        starts[cnt] = vertices[1, a]
                        stops[cnt] = vertices[0, b]
                        cnt += 1
        for m in range(n):
            if not excluded[m] and vertices[1, m] > vertices[0, m]:
                starts[cnt] = vertices[0, m]
                stops[cnt] = vertices[1, m]
                cnt += 1
        ### after sorting by start, two intervals overlap, if and only if one starts before its predecessor ends
        order = np.argsort(starts[:cnt])
        for q in range(cnt - 1):
            if starts[order[q + 1]] < stops[order[q]]:
                overlapping = True
                break

    ### if at any position more than one exon or intron -> is_alt_spliced
    if overlapping:
        return (0, 1, 1)
    ### if not alt_spliced but term_alt or init_alt --> is_alt
    if n_alt == 0:
        return (0, 0, 0)
    return (0, 0, 1)


class Gene:
    
    def __init__(self, name=None, start=None, stop=None, chr=None, strand=None, source=None, gene_type=None, gene_symbol=None):
        self.name = name
        self.start = start
        self.stop = stop
        self.exons = []
        self.introns_anno = set([])
        self.chr = chr
        if strand in ['+', '-']:
            self.strand = strand
        else:
            warnings.warn('WARNING: strand of gene was provided as %s - automatically set to +' % str(strand))
            self.strand = '+'
        self.transcripts = []
        self.source = source
        self.splicegraph = Splicegraph()
        self.segmentgraph = Segmentgraph()
        self.gene_type = gene_type
        self.is_alt = None
        self.is_alt_spliced = None
        self.symbol = gene_symbol

    def __eq__(self, other):
        return isinstance(other, Gene) and \
               (self.name == other.name) and \
               (self.start == other.start) and \
               (self.stop == other.stop) and \
               (len(self.exons) == len(other.exons)) and \
               min([np.all(self.exons[i] == other.exons[i]) for i in range(len(self.exons))]) and \
               (self.chr == other.chr) and \
               (self.strand == other.strand) and \
               (self.source == other.source) and \
               (self.transcripts == other.transcripts) and \
               (self.splicegraph == other.splicegraph) and \
               (self.segmentgraph == other.segmentgraph) and \
               (self.gene_type == other.gene_type) and \
               (self.is_alt == other.is_alt) and \
               (self.is_alt_spliced == other.is_alt_spliced) and \
               (self.symbol == other.symbol)


    def add_exon(self, exon, idx):
        if idx > (len(self.exons) - 1): 
            self.exons.append(np.zeros((0, 2), dtype='int'))
        self.exons[idx] = np.concatenate((self.exons[idx], np.asarray(exon)[np.newaxis, :]))

    def label_alt(self):
        """Labels the gene as alternatively spliced (is_alt_spliced) and as having alternative transcript starts or ends (is_alt)."""

        vertices = self.splicegraph.vertices
        edges = self.splicegraph.edges
        ### the compiled version takes the usual graphs, everything else is left to the general implementation
        if (type(vertices) is np.ndarray and type(edges) is np.ndarray and vertices.dtype == np.int64 and edges.dtype in (np.int64, np.bool_) and
                vertices.ndim == 2 and vertices.shape[0] == 2 and edges.ndim == 2 and edges.shape == (vertices.shape[1], vertices.shape[1]) and
                2 <= vertices.shape[1] <= _MAX_COMPILED_EXONS):
            (status, is_alt_spliced, is_alt) = _label_alt_kernel(vertices, edges)
            if status == 0:
                self.is_alt_spliced = is_alt_spliced
                self.is_alt = is_alt
                return
        self._label_alt_general()

    def _label_alt_general(self):
        """label_alt for any splice graph; graphs with exons that end before they start are handled like the original implementation does it."""

        vertices = self.splicegraph.vertices
        if type(vertices) is np.ndarray and vertices.ndim == 2 and vertices.shape[0] == 2 and not np.any(vertices[1, :] < vertices[0, :]):
            self._label_alt_numpy()
        else:
            self._label_alt_original()

    def _label_alt_numpy(self):

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
        
        ### only the upper triangle of the adjacency matrix is looked at
        upper = np.triu(edges != 0, 1)
        has_in = upper.any(axis=0)
        has_out = upper.any(axis=1)
        ### exons that share their end / start with another exon
        _, inverse, counts = np.unique(vertices[1, :], return_inverse=True, return_counts=True)
        shares_end = counts[inverse.ravel()] > 1
        _, inverse, counts = np.unique(vertices[0, :], return_inverse=True, return_counts=True)
        shares_start = counts[inverse.ravel()] > 1

        init_alt_idx = []
        term_alt_idx = []
        ### all other exons neither are initial nor terminal alternative exons
        for i in np.flatnonzero((~has_in & shares_end) | (~has_out & shares_start)):
            i = int(i)
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
        excluded = np.zeros((num_exons,), dtype='bool')
        excluded[np.array(init_alt_idx + term_alt_idx, dtype='int')] = True
        take_idx = np.where(~excluded)[0]
        if take_idx.shape[0] > 0:
            vertices = self.splicegraph.vertices[:, take_idx]
            edges = self.splicegraph.edges[take_idx, :][:, take_idx]
          
            start = vertices[0, :].min()
            stop = vertices[1, :].max()
            if stop <= start:
                ### as before, an empty array has no maximum
                np.zeros((stop - start,)).max()

            ### Is any position covered by more than one intron or exon?  Introns are the stretches between the end of an
            ### exon and the start of a following exon it is connected to.  Instead of counting the coverage of every position
            ### of the gene, look for two (half open) intervals that overlap.  After sorting by start this is the case, if
            ### and only if an interval starts before its predecessor ends: if an interval overlaps any earlier one, the
            ### interval right after that earlier one starts inside it as well.
            edge_from, edge_to = np.nonzero(np.triu(edges == 1, 1))
            interval_start = np.concatenate((vertices[1, edge_from], vertices[0, :]))
            interval_stop = np.concatenate((vertices[0, edge_to], vertices[1, :]))
            non_empty = interval_stop > interval_start
            interval_start = interval_start[non_empty]
            interval_stop = interval_stop[non_empty]
            order = np.argsort(interval_start, kind='stable')
            interval_start = interval_start[order]
            interval_stop = interval_stop[order]
            overlapping = bool(np.any(interval_start[1:] < interval_stop[:-1]))
        else:
            overlapping = False
      
        ### if at any position more than one exon or intron -> is_alt_spliced
        if overlapping:
            self.is_alt_spliced = 1
            self.is_alt = 1
        else:
            self.is_alt_spliced = 0
            ### if not alt_spliced but term_alt or init_alt --> is_alt
            if len(init_alt_idx) == 0 and len(term_alt_idx) == 0:
                self.is_alt = 0 
            else:
                self.is_alt = 1 

    def _label_alt_original(self):
        """The original implementation, which counts the coverage of every position of the gene (slow; the reference of the others)."""

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

    def get_non_alt_seg_ids(self):
        """Returns a list of segment indices that are not alternatively used in the graph.
           
           Currently, the definition of non-alternative means the segment is not spanned
           by any intron. In the narrow sense of non-alternative, one should also exclude
           alternative starts and ends but this could lead to cases where no segment is 
           considered as non-alternative, so we use this approximation for now.
        """

        tmp = np.ones((self.segmentgraph.seg_edges.shape[0],), dtype='bool')
        for i in range(self.segmentgraph.seg_edges.shape[0] - 1):
            ### get index of last acceptor
            idx = np.where(self.segmentgraph.seg_edges[i, i + 1:])[0]
            ### mask all segments between current segment and acceptor
            if idx.shape[0] > 0:
                tmp[i + 1:idx[-1] + i + 1] = 0

        return np.where(tmp)[0]

    def to_sparse(self):
        
        self.splicegraph_edges_shape = self.splicegraph.edges.shape
        (self.splicegraph_edges_data, self.splicegraph_edges_row, self.splicegraph_edges_col) = _matrix_to_coo(self.splicegraph.edges)
        self.splicegraph.edges = None

        if hasattr(self, 'edge_count'):
            self.edge_count_shape = self.edge_count.shape
            (self.edge_count_data, self.edge_count_row, self.edge_count_col) = _matrix_to_coo(self.edge_count)
            self.edge_count = None 
        

    def from_sparse(self):
        
        if hasattr(self, 'splicegraph_edges_data'):
            self.splicegraph.edges = _matrix_from_coo(self.splicegraph_edges_data, self.splicegraph_edges_row, self.splicegraph_edges_col, self.splicegraph_edges_shape)
            del self.splicegraph_edges_data
            del self.splicegraph_edges_shape
            del self.splicegraph_edges_row
            del self.splicegraph_edges_col

        if hasattr(self, 'edge_count_data'):
            self.edge_count = _matrix_from_coo(self.edge_count_data, self.edge_count_row, self.edge_count_col, self.edge_count_shape)
            del self.edge_count_data
            del self.edge_count_shape
            del self.edge_count_row
            del self.edge_count_col


    def populate_annotated_introns(self):
        self.introns_anno = set()
        for tx in self.exons:
            sidx = np.argsort([_[0] for _ in tx])
            for e in range(1, len(sidx)):
                self.introns_anno.add((tx[sidx[e - 1]][1], tx[sidx[e]][0]))
                
