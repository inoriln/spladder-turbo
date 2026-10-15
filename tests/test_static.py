"""Static checks of the package."""

import builtins
import dis
import importlib
import pkgutil
import types

import spladder

### unfinished code of the original (it cannot run: the names it uses do not exist anywhere), nothing calls it with the options available
KNOWN_UNFINISHED = {
    ('spladder.alt_splice.write', 'write_events_txt', 'load'),                    # annotation of the events (anno_fn): never ported
    ('spladder.alt_splice.write', 'write_events_txt', 'sort'),
    ('spladder.alt_splice.write', 'write_events_txt', 'genes'),
    ('spladder.core.gen_graphs', 'gen_graphs', 'infer_splice_graph'),             # --re-infer-sg: the function does not exist
    ('spladder.alt_splice.quantify', 'quantify_mutex_exons', 'cov'),               # the per-event quantification (quantify_from_graph) is not called by anything and does not run with Python 3
    ('spladder.editgraph', 'infer_splice_graph_caller', 'test_exon_id'),
    ('spladder.editgraph', 'infer_splice_graph_caller', 'Tru'),
    ('spladder.reads', 'get_all_data_uncollapsed', 'lenfilenames'),
    ('spladder.reads', 'get_all_data_uncollapsed', 'ot'),
}


def _code_objects(code):
    yield code
    for const in code.co_consts:
        if isinstance(const, types.CodeType):
            yield from _code_objects(const)


def _undefined_global_names():
    """(module, function, name) of the names a function loads as globals that its module does not define (star imports included)."""

    found = set()
    for info in pkgutil.walk_packages(spladder.__path__, 'spladder.'):
        try:
            module = importlib.import_module(info.name)
        except ImportError:
            continue
        path = getattr(module, '__file__', None)
        if not path or not path.endswith('.py'):
            continue
        with open(path) as fd:
            top = compile(fd.read(), path, 'exec')
        codes = [c for c in _code_objects(top) if c is not top]
        ### (names that a function assigns with `global`)
        assigned = {ins.argval for c in codes for ins in dis.get_instructions(c) if ins.opname == 'STORE_GLOBAL'}
        for code in codes:
            for ins in dis.get_instructions(code):
                if ins.opname in ('LOAD_GLOBAL', 'LOAD_NAME') and ins.argval not in vars(module) and ins.argval not in assigned and not hasattr(builtins, ins.argval):
                    found.add((info.name, code.co_name, ins.argval))
    return found


def test_the_functions_of_the_package_use_names_that_exist():
    undefined = _undefined_global_names()
    assert undefined - KNOWN_UNFINISHED == set()
    ### (the list does not outlive the code it excuses)
    assert KNOWN_UNFINISHED - undefined == set()
