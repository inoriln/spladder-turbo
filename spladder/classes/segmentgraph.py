import numpy as np

class Segmentgraph:

    def __init__(self, gene = None):
    
        self.segments = np.zeros((2, 0), dtype='int')
        self.seg_match = np.zeros((0, 0), dtype='bool')
        self.seg_edges = np.zeros((0, 0), dtype='bool')

        if gene is not None:
            self.from_gene(gene)

    def __eq__(self, other):
        return isinstance(other, Segmentgraph) and \
               np.all(self.segments.shape == other.segments.shape) and \
               np.all(self.segments == other.segments) and \
               np.all(self.seg_match.shape == other.seg_match.shape) and \
               np.all(self.seg_match == other.seg_match) and \
               np.all(self.seg_edges.shape == other.seg_edges.shape) and \
               np.all(self.seg_edges == other.seg_edges)

    def is_empty(self):
        
        return (self.segments.shape[1] == 0) and (self.seg_match.shape[0] == 0) and (self.seg_edges.shape[0] == 0)

    def from_gene(self, gene): 

        sg = gene.splicegraph.vertices
        breakpoints = np.unique(sg.ravel())

        ### a stretch between two breakpoints is a segment, if more exons started than ended before its end
        n_started = np.searchsorted(np.sort(sg[0, :]), breakpoints[1:], side='left')
        n_ended = np.searchsorted(np.sort(sg[1, :]), breakpoints[1:], side='left')
        is_segment = n_started > n_ended
        self.segments = np.vstack([breakpoints[:-1][is_segment], breakpoints[1:][is_segment]]).astype('int')

        ### match nodes to segments
        self.seg_match = ((sg[0, :][:, np.newaxis] <= self.segments[0, :][np.newaxis, :]) & (sg[1, :][:, np.newaxis] >= self.segments[1, :][np.newaxis, :]))

        ### create edge graph between segments
        self.seg_edges = np.zeros((self.segments.shape[1], self.segments.shape[1]), dtype='bool')
        k, l = np.where(np.triu(gene.splicegraph.edges))

        has_segment = self.seg_match.any(axis=1)
        if k.shape[0] == 0:
            pass
        elif has_segment[k].all() and has_segment[l].all():
            ### donor segment: last segment of the donor exon, acceptor segment: first segment of the acceptor exon
            d = self.seg_match.shape[1] - 1 - np.argmax(self.seg_match[k, ::-1], axis=1)
            a = np.argmax(self.seg_match[l, :], axis=1)
            self.seg_edges[d, a] = True
        else:
            ### an exon without any segment: the loop below raises the same error as before
            for m in range(k.shape[0]):
                d = np.where(self.seg_match[k[m], :])[0][-1]
                a = np.where(self.seg_match[l[m], :])[0][0]
                self.seg_edges[d, a] = True
