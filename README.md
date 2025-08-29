# Towards provable energy stable overset grid methods using sub-cell summation-by-parts operators

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/TODO/zenodo.TODO.svg)](https://doi.org/TODO/zenodo.TODO)


This repository contains information and code to reproduce the results presented in the
article
```bibtex
@online{glaubitz2025towards,
  title={Towards provable energy stable overset grid methods using
         sub-cell summation-by-parts operators},
  author={Glaubitz, Jan and Lampert, Joshua and Winters, Andrew R and Nordström, Jan},
  year={2025},
  month={09},
  eprint={TODO},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{glaubitz2025towardsRepro,
  title={Reproducibility repository for
         "{T}owards provable energy stable overset grid methods using
         sub-cell summation-by-parts operators"},
  author={Glaubitz, Jan and Lampert, Joshua and Winters, Andrew R and Nordström, Jan},
  year={2025},
  howpublished={\url{https://github.com/JoshuaLampert/2025\_overset\_grid\_sub-cell}},
  doi={TODO/zenodo.TODO}
}
```

## Abstract

Overset grid methods handle complex geometries by overlapping simpler, geometry-fitted grids to cover the
original, more complex domain. However, ensuring the stability of the procedures used to couple these
grids—particularly at high orders—remains a practical and theoretical challenge. In this work, we address
this gap by developing a discrete counterpart to the recent well-posedness analysis of Kopriva, Gassner, and
Nordström for continuous overset domain initial-boundary-value problems (IBVPs). To this end, we intro-
duce the novel concept of sub-cell summation-by-parts (SBP) operators. These discrete derivative operators
mimic integration by parts at a sub-cell level. By exploiting this sub-cell SBP property, we develop provably
conservative and energy-stable overset grid methods for stationary meshes, thereby resolving longstanding
stability issues in the field.


## Numerical experiments

To reproduce the numerical experiments presented in this article, you need
to install [Julia](https://julialang.org/). The numerical experiments presented
in this article were performed using Julia v1.11.6.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface. Then, you need to start
Julia in the `code` directory of this repository and follow the instructions
described in the `README.md` file therein.


## Authors

- Jan Glaubitz (Linköping University, Sweden)
- Joshua Lampert (University of Hamburg, Germany)
- Andrew R. Winters (Linköping University, Sweden)
- Jan Nordström (Linköping University, Sweden and University of Johannesburg, South Africa)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
