#! /usr/bin/env python
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2021 Andre Kahles, Jonas Behr, Gunnar Raetsch
# Copyright (C) 2009-2011 Max Planck Society
# Copyright (C) 2012-2016 Memorial Sloan-Kettering Cancer Center
# Copyright (C) 2016-2021 ETH Zurich

import copy
import sys
import os
import numpy as np
import pickle
import h5py

from . import parallel
from . import settings
from .core.spladdercore import spladder_core
from .alt_splice.collect import collect_events
from .alt_splice.analyze import analyze_events
from .alt_splice.verify import graph_filenames
from .count import count_graph_coverage_wrapper, collect_single_quantification_results, compute_gene_expression
from .editgraph import filter_by_edgecount
from .merge import run_merge
from .helpers import *
from .utils import dump_pickle, load_pickle, load_pickle_shared

from .spladder_prep import prep_annotation, prep_sparse_bam_filtered, prep_sparse_bam_full

from .classes.gene import Gene
sys.modules['modules.classes.gene'] = Gene

### the event types collected by one job each, jobs of about the same duration; the two alternative prime ends come from the same detection
_EVENT_GROUPS = [['exon_skip'], ['mult_exon_skip'], ['intron_retention', 'alt_3prime', 'alt_5prime', 'mutex_exons']]


def _event_collection_jobs(options):
    """The options for the jobs that collect the events (one job for unusual event types)."""

    known = [t for group in _EVENT_GROUPS for t in group]
    if any(t not in known for t in options.event_types):
        return [options]
    jobs = []
    for group in _EVENT_GROUPS:
        selected = [t for t in group if t in options.event_types]
        if len(selected) > 0:
            job_options = copy.copy(options)
            job_options.event_types = selected
            jobs.append(job_options)
    return jobs


def _analyze_task(number):
    """Verifies and reports the events of one type (the work of a process of _analyze_events)."""

    tasks, options = parallel.shared()
    idx, event_type = tasks[number]
    if options.merge == 'single':
        analyze_events(event_type, options.bam_fnames, options, sample_idx=idx)
    else:
        analyze_events(event_type, options.bam_fnames, options)


def _analyze_events(idxs, options):
    """Verifies and reports the events of all types; the types are independent of each other and are done at the same time if there are processes to do so."""

    tasks = [(idx, event_type) for idx in idxs for event_type in options.event_types]
    ### (one splice graph is shared by the tasks: not for the graphs of many samples)
    if len(idxs) != 1 or len(tasks) < 2 or not parallel.can_run_parallel(options, len(tasks)):
        for (idx, event_type) in tasks:
            if options.merge == 'single':
                analyze_events(event_type, options.bam_fnames, options, sample_idx=idx)
            else:
                analyze_events(event_type, options.bam_fnames, options)
        return

    ### what the tasks share: the splice graph (read only)
    ### (the index of the counts is not made here: it is only made by a task that has events to verify, as in a serial run)
    if options.verify_alt_events:
        for idx in idxs:
            fn_genes = graph_filenames(options, idx)[0]
            if os.path.exists(fn_genes):
                try:
                    load_pickle_shared(fn_genes)
                except Exception:
                    ### (a graph that cannot be read makes the tasks that need it fail, as in a serial run)
                    pass
    for _ in parallel.run_shards(_analyze_task, (tasks, options), len(tasks), len(tasks), daemon=False):
        pass


def _prep_workdir(options):

    ### create out-directory, spladder and tmp sub-directory
    ### (exist_ok: several SplAdder processes may be started on the same, still empty output directory)
    os.makedirs(options.outdir, exist_ok=True)
    os.makedirs(os.path.join(options.outdir, 'spladder'), exist_ok=True)
    os.makedirs(options.tmpdir, exist_ok=True)


def spladder(options):

    with parallel.sigterm_unwinds():
        try:
            return _spladder(options)
        finally:
            ### a job that is still running is of no use any more, e.g., after an error
            parallel.abort_background()


def _spladder(options):

    ### parse parameters from options object
    options = settings.parse_args(options)

    ### load confidence level settings
    if not options.no_reset_conf:
        options = settings.set_confidence_level(options)

    ### do not compute components of merged set, if result file already exists
    fn_out_merge = get_filename('fn_out_merge', options)
    fn_out_merge_val = get_filename('fn_out_merge_val', options)

    ### preparatory work
    _prep_workdir(options)
    options = prep_annotation(options)

    ###########################################################################
    ### GRAPH GENERATION
    ###########################################################################

    if not os.path.exists(fn_out_merge):

        ### convert input BAMs to sparse arrays - filtered case
        if options.sparse_bam:
            prep_sparse_bam_filtered(options)

        ### handle single graphs
        if options.merge in ['single', 'merge_graphs', 'merge_all']:
            graph_fname = lambda idx: '%s/spladder/genes_graph_conf%i.%s.pickle' % (options.outdir, options.confidence, options.samples[idx])
            ### samples without a graph (a sample name occurring twice is computed once, as in the serial loop)
            todo = []
            for idx in range(len(options.samples)):
                if not os.path.exists(graph_fname(idx)) and not [i for i in todo if graph_fname(i) == graph_fname(idx)]:
                    todo.append(idx)
            if parallel.can_run_parallel(options, len(todo)):
                parallel.build_sample_graphs(options, todo, graph_fname, spladder_core)
            else:
                for idx in range(len(options.samples)):
                    out_fname = graph_fname(idx)
                    if not os.path.exists(out_fname):
                        spladder_core(options.bam_fnames[idx], out_fname, options)
        ### handle merged bams
        if options.merge in ['merge_bams', 'merge_all']:
            out_fname = '%s/spladder/genes_graph_conf%i.merge_bams.pickle' % (options.outdir, options.confidence)
            if not os.path.exists(out_fname):
                spladder_core(options.bam_fnames, out_fname, options)
        ### handle merged graphs
        if options.merge in ['merge_graphs', 'merge_all']:
            run_merge(options.samples, options)

    ### generate validated version of splice graph
    if options.merge == 'merge_graphs' and options.validate_sg and not os.path.exists(fn_out_merge_val):
        (genes, inserted) = load_pickle(fn_out_merge)
        genes = filter_by_edgecount(genes, options)
        dump_pickle((genes, inserted), fn_out_merge_val)
        del genes

    ###########################################################################
    ### GRAPH QUANTIFICATION
    ###########################################################################

    ### convert input BAMs to sparse arrays - unfiltered case
    if options.sparse_bam:
        prep_sparse_bam_full(options)

    if options.merge == 'single' or options.qmode == 'collect':
        idxs = list(range(len(options.samples)))
    else:
        idxs = [0]

    ### collecting the events needs the splice graph only (not the counts): it is done while the graph is quantified
    events_jobs = None
    if options.extract_as and options.quantify_graph and parallel.can_run_parallel(options, 2):
        events_jobs = [parallel.Background(collect_events, job_options) for job_options in _event_collection_jobs(options)]

    for idx in idxs:
        ### don't do any counting in chunked merge mode, unless we are in the final chunk level
        if len(options.chunked_merge) > 0:
            curr_level, max_level, chunk_start, chunk_end = [int(_) for _ in options.chunked_merge[0]]
            if curr_level < max_level:
                break

        ### get count output file
        if options.merge == 'single':
            fn_in_count = get_filename('fn_count_in', options, options.samples[idx])
            fn_out_count = get_filename('fn_count_out', options, options.samples[idx])
        elif options.merge == 'merge_graphs' and options.qmode == 'single':
            fn_in_count = get_filename('fn_count_in', options)
            fn_out_count = get_filename('fn_count_out', options, options.samples[0])
        else:
            fn_in_count = get_filename('fn_count_in', options)
            fn_out_count = get_filename('fn_count_out', options)
        fn_out_gene_count = re.sub(r'.count.hdf5$', '', fn_out_count) + '.gene_exp.hdf5'

        ### count segment graph
        if options.quantify_graph:
            if not os.path.exists(fn_out_count):
                ### collect graph counts
                if options.merge == 'single':
                    count_graph_coverage_wrapper(fn_in_count, fn_out_count, options.bam_fnames, options, sample_idx=idx)
                elif options.merge == 'merge_graphs' and options.qmode == 'single':
                    count_graph_coverage_wrapper(fn_in_count, fn_out_count, options.bam_fnames, options, qmode='single')
                elif options.merge == 'merge_graphs' and options.qmode == 'collect':
                    collect_single_quantification_results(fn_out_count, idxs, options)
                else:
                    count_graph_coverage_wrapper(fn_in_count, fn_out_count, options.bam_fnames, options)
                
            ### collect gene expression counts from graph
            if not os.path.exists(fn_out_gene_count):
                if options.merge == 'single':
                    compute_gene_expression(options, fn_in_count, fn_out_count, fn_out_gene_count, sample_idx=[0])
                else:
                    compute_gene_expression(options, fn_in_count, fn_out_count, fn_out_gene_count)

    ###########################################################################
    ### ALT EVENT COLLECTION
    ###########################################################################

    ### handle alternative splicing part
    if options.extract_as:
        if events_jobs is None:
            collect_events(options)
        else:
            for job in events_jobs:
                job.join()

        if options.quantify_graph:
            _analyze_events(idxs, options)


