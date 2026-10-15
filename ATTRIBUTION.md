# Attribution

**spladder-turbo** is an unofficial, performance-optimized fork of
[SplAdder](https://github.com/ratschlab/spladder) 3.1.1 (upstream commit `65ceec8`). It is not
affiliated with, endorsed by, or supported by the SplAdder authors or their institutions.

The method, the definitions of the splicing events and the output formats are SplAdder's. This fork
changes how the work is done (reading the alignments, scheduling of parallel processes, compiled
kernels) and fixes a number of defects; see [docs/PERFORMANCE_CHANGES.md](docs/PERFORMANCE_CHANGES.md).
The git history of this repository contains the complete history of SplAdder.

## Original authors and copyright holders

SplAdder was written by Cheng Soon Ong, Gunnar Raetsch and Andre Kahles (see [AUTHORS](AUTHORS)),
with contributions from others (see the git history). All copyright notices of SplAdder are retained
unchanged: [COPYRIGHT](COPYRIGHT) and the headers of the source files, which also name the Max Planck
Society, the Memorial Sloan-Kettering Cancer Center and ETH Zurich.

## License

- [COPYRIGHT](COPYRIGHT) (BSD 3-Clause) is the license under which SplAdder is distributed upstream.
  [LICENSE](LICENSE) contains the same text, with the notice for the modifications of this fork appended.
  As the third clause says, the names of the copyright holders and contributors are not used here to
  endorse or promote this fork.
- `spladder/spladder_build.py` additionally carries an upstream notice of the GNU General Public License,
  version 3 or later. It is retained unchanged. The license information of upstream is not fully
  consistent on this point (COPYRIGHT and `setup.py` say BSD).
- Files that were added by this fork are distributed under the BSD 3-Clause license as well.

## Citation

If you use this software, please cite the SplAdder publication:

Kahles A, Ong CS, Zhong Y, Rätsch G. SplAdder: identification, quantification and testing of alternative
splicing events from RNA-Seq data. *Bioinformatics* 2016;32(12):1840-1847.
doi:[10.1093/bioinformatics/btw076](https://doi.org/10.1093/bioinformatics/btw076)
