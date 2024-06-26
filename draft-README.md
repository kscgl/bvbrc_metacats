# Metadata-driven Comparative Analysis Tool (meta-CATS)

## Overview

The meta-CATS metadata genome comparison tool takes sequence data and determines the aligned positions that significantly differ between two (or more) user-specified groups. Once an analysis is started, a multiple sequence alignment is performed if the input was unaligned (such as from a database query). A chi-square test of independence is then performed on each non-conserved column of the alignment, to identify those that have a non-random distribution of bases. A quantitative statistical value of variation is computed for all positions. Columns that are perfectly conserved will not be identified as statistically significant. All other non-conserved columns will be evaluated to determine whether the p-value is lower than the specified threshold value. Terminal gaps flanking the aligned sequences will not be taken into account for the analysis.  



## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

This module provides the following application specfication(s):
* [MetaCATS](app_specs/MetaCATS.md)


## See also

* [Metadata-driven Comparative Analysis Tool (meta-CATS) Quick Reference](https://www.bv-brc.org/docs/quick_references/services/metacats.html)
* [Metadata-driven Comparative Analysis Tool (https://www.bv-brc.org/docs/meta-CATS.html)](https://bv-brc.org/app/MetaCATS)
* [Metadata-driven Comparative Analysis Tool (https://www.bv-brc.org/docs/meta-CATS.html) Tutorial](/tutorial/metacats/metacats)



## References

