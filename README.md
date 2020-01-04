# BioGraph 

![](logo.png)

BioGraph is a library for manipulating structures in PDB files as graphs,
with the purpose of training machine learning models. Additionally
it provides functionality for making pipelines based on those graphs.
It is based on [BioPython](https://github.com/biopython/biopython).

In short, this package provides the following functionality:

* representing a PDB file as a `networkx.Graph` through a structure
generator. Different generators can be used to create different
graphs. For example, one can generate a static contacts graph,
or a graph that captures the topology of the protein (e.g.
Delauney filtration).

* adding features to either the protein (represented as a pandas
dataframe) or the graph through alignment.

* diffusing or propagating features along neighbors and optionally
converting the end result to a dataframe for use in machine learning
models that a priori take the graph structure as input.

* making cross validation fold according to CDHit groups to avoid
bias in training.


### Installation

Since BioGraph has quite a few dependencies that need to be compiled
prior to installing, we've yet to figure out if it's possible only
using `pip`.

Fortunately `setuptools` provides more than enough functionality
to automate these tasks, so just running:

```
# python setup.py install
```

should do the trick.

### Examples

Example notebooks are provided in `examples/`. The most
straightforward machine learning example is in
`examples/ExamplePipeline.ipynb`.
