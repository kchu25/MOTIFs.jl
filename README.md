# MOTIFs

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/MOTIFs.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/MOTIFs.jl/dev/)
[![Build Status](https://github.com/kchu25/MOTIFs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/MOTIFs.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/MOTIFs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/MOTIFs.jl)

General purpose motif discovery package that includes the discovery of flexible (long or gapped) motifs.

This motif discovery package is currently being actively developed and the method is currently under review with a journal.

# Installation
To install MOTIFs.jl use Julia's package manager:
```
pkg> add MOTIFs
```

# Usage
````julia
using MOTIFs

# Do motif discovery on a set of DNA sequences in a fasta file, where the `<fasta-path>` and `<output-folder-path>` are the absolute filepaths as strings.

discover_motifs(<fasta-path>, <output-folder-path>)

````

# Software requirements 
 This package currectly requires [Weblogo](http://weblogo.threeplusone.com/manual.html#download). You will need python3 and install Weblogo with following command:
 ```bash
 pip3 install weblogo
 ```
