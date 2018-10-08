PHY982 project
==============

**Quick links**: [paper][a], [presentation][p] ([pdf][f]).

Source code for my [PHY982][1] class project.

Building
--------

First off, make sure you have the necessary [dependencies](#dependencies).

  - To build the paper and presentation:

        make paper.pdf index.html

  - To run the calculations:

        make run

Dependencies
------------

### Figures (required for both paper and presentation)

  - [Inkscape](https://inkscape.org) (tested on 0.91)

### Paper

  - [LaTeX](http://latex-project.org)
    with the following packages:

      - bm
      - graphicx
      - hyperref
      - revtex4-1
      - xcolor

### Presentation

  - [Pandoc](https://pandoc.org) (tested on 1.13.2.1)

### Calculations

  - [Lapack](https://netlib.org/lapack)

  - [Python](https://python.org) 3.4+ or 2.7+

      - [matplotlib](https://matplotlib.org)
      - [numpy](https://numpy.org)
      - [scipy](https://scipy.org)

  - [Slatec](https://netlib.org/slatec) (tested on 4.1)

[1]: https://people.nscl.msu.edu/~nunes/phy982/phy982web2015.htm
[a]: https://xrf.github.io/phy982-proj/paper.pdf
[p]: https://xrf.github.io/phy982-proj
[f]: https://xrf.github.io/phy982-proj/proj.pdf
