## White et al. 2012 Analysis

The code and data in this repository allow the analyses in White et al. 2012
([http://dx.doi.org/10.1890/11-2177.1](http://dx.doi.org/10.1890/11-2177.1)) to
be fully replicated for four of the six datasets Breeding Bird Survey, Mammal
Community Database, Forest Inventory and Analysis, and Alwyn Gentry's Tree
Transects). The other two datasets were obtained under agreements restricting
the publication of raw data, but simulation results and figures can still be
generated for these datasets.

### Setup

Requirements: Python 2.x and the following Python modules: numpy, scipy,
matplotlib, mpmath, mpl_toolkits (for figures), and mpl_toolkits.basemap (for
figures). You will also need two of our custom Python modules: METE (at the
white-etal-2012 tag; https://github.com/weecology/METE) and macroecotools (at
the white-etal-2012 tag; https://github.com/weecology/macroecotools). These
modules can be installed by running the following commands from the command line
(with sufficient permissions):

```sh
git clone https://github.com/weecology/METE.git
cd METE
git checkout white-etal-2012
python setup.py install
cd ..
git clone https://github.com/weecology/macroecotools.git
cd macroecotools
git checkout white-etal-2012
python setup.py install
```

### Replicate analyses

The analyses can be replicated by running the following commands from the
command line.

Run all analyses and generate figures: `python mete_sads.py ./data/ all 100`
(where 100 is the number of simulations you wish to conduct and can be replaced
with any positive integer)


Run portions of the analysis pipeline:

* Empirical analyses: `python mete_sads.py ./data/ empir`
* Simulation analyses: `python mete_sads.py ./data/ sim 100` (where 100 is the
  number of simulations you wish to conduct and can be replaced with any
  positive integer)
* Figures: `python mete_sads.py ./data/ figs`

On Windows `./data/` should be replaced with `.\data\` to match the relevant
path conventions.

The data included in the repository shows the data, including intermediate
steps, that were used in the paper, and can be used to reproduce pieces of the
pipeline. So, if you just want to remake the figures, the code will use the
provided data. However, we also want to be able to rerun the entire
pipeline. So, if you start at the beginning and run with the 'all' argument, the
code will start with just the rawest form of the data and recreate every
intermediate step including recreating the intermediate data files, exactly as
was done to create them for the paper in the first place. This will delete the
existing data files and replace them with the recalculated versions. Due to some
minor post-publication bug fixes the exact form of the resulting data files will
differ slightly from those currently in the repository. None of these
differences influence the results of the paper.

Please note that these analyses involve both a large amount of data and a lot of
computational work and therefore take a long time to run. Expect the empirical
analysis to take up to a day, and simulations to take up to a week on an 8-core
server. Generating figures takes about one hour due to the neighborhood
calculations required for the color ramps on the observed-predicted plots.

### Data use

Data is provided in this supplement for the purposes of replication and is not
presented in such a way as to be generally useful for additional analyses. If
you wish to use these datasets for additional research they should be obtained
from relevant data providers. For BBS, MCDB, FIA, and Gentry this can be done
automatically by using the EcoData Retriever
(http://ecodataretriever.org).

The ``mete_sads_data.py`` script is intendend to allow for downloading the raw data
from the original source, importing it into MySQL and then executing the queries
we used to get the data in the ./data folder. However, it is not currently
complete. Hopefully, we'll get back to it soon. Sorry for the inconvenience.

### License (MIT)

Copyright (c) 2012 Weecology

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
