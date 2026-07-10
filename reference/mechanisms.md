# Mechanism base class and the Gaussian-process implementation

A `Mechanism` encapsulates the probabilistic model for a single node's
conditional distribution `P(X | parents)`. Kagu's sole built-in
mechanism is the Gaussian process
([GPMechanism](https://causalabs.github.io/kagu-r/reference/GPMechanism.md)),
a flexible nonparametric default. To add a new backend, subclass
`Mechanism` and implement the methods below.
