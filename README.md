# Finding Motifs Using DNA Images Derived From Sparse Representations

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/MOTIFs.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/MOTIFs.jl/dev/)
[![Build Status](https://github.com/kchu25/MOTIFs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/MOTIFs.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/MOTIFs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/MOTIFs.jl)

General purpose motif discovery package that includes the discovery of flexible (long or gapped) motifs.

# Motivation
(coming soon)

# Installation
To install MOTIFs.jl use Julia's package manager:
```
pkg> add MOTIFs
```

# Usage
````julia
using MOTIFs

# Do motif discovery on a set of DNA sequences in a fasta file, 
# where the `<fasta-path>` and `<output-folder-path>` are the 
# absolute filepaths as strings.

discover_motifs(<fasta-path>, <output-folder-path>)

# for example

discover_motifs("home/shane/mydata/fasta.fa", 
                "home/shane/mydata/out/")
````

# Interpret the results
(coming soon)


# Software requirements 
 This package currectly requires [Weblogo](http://weblogo.threeplusone.com/manual.html#download) for PWM plotting. Install Weblogo by running the following command with python3 and pip3:
 ```bash
 pip3 install weblogo
 ```

# Hardware requirements
For now, a GPU is required for this package due to the use of [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) to accelerate some computations. I plan to implement a CPU extension in the future.


# Adjustable Hyperparameters
````julia

# The user can adjust the number of epochs for training the network.
discover_motifs(<fasta-path>, <output-folder-path>; num_epochs=10)

````

# Citation <a name="cite"></a>

The paper presenting this method has been published in [Oxford Bioinformatics](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btad378/7192989?utm_source=advanceaccess&utm_campaign=bioinformatics&utm_medium=email). It can be cited using the following BibTex entry:
```
@article{chu2023finding,
  title={Finding Motifs Using DNA Images Derived From Sparse Representations},
  author={Chu, Shane K and Stormo, Gary D},
  journal={Bioinformatics},
  pages={btad378},
  year={2023},
  publisher={Oxford University Press}
}
```

# Contact

If you have any questions or suggestions regarding the usage or source code, please feel free to reach out to me at <skchu@wustl.edu>.