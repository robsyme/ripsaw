# Ripsaw #

This is a tool that is (at the moment) just OcculterCut with some bug
fixes.  Perhaps one day it will grow into something that can stand on
its own, who knows?


## Ripsaw vs OcculterCut

The largest 'bugfix' is that OcculterCut predicts GC and AT regions
over scaffold gaps, which may skew results of an inattentive
user. This tool also parallelises the calculations over the available
cores.

More deeply, Ripsaw aims to expand the information it uses for
sequence classification. OcculterCut calculates the Jensen–Shannon
divergence of a binomical distribution (A+T vs G+C). At present,
Ripsaw calculates the J-S divergence of the 14-value distribution of
all possible dinucleotides. It finds breakpoints that maximise the
dinucleotide distribution divergence (in the same manner as
OcculterCut) and then outputs the small bins with their 14-value
distibution calculations. From there, I need to classify the bins and
then re-join neighboring bins that have the same class.

## Classification

OcculterCut fits a mixture of Cauchy Distributions to the GC density
distribution and finds a value that most cleanly intersects the
distributions. There are other possibilities (k-means, SOM, etc) which
I have only made tentative attempts at here.

### How do I get set up? ###

Installation/compilation is 
- Clone the repository
- `go build .`

### Contribution guidelines ###

* Writing tests
Yeah, I'll get to that eventually. At the moment the project is too simple, 
but that might change.

### Who do I talk to? ###

Any questions - talk to me - rob.syme@curtin.edu.au
