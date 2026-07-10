# kagu: Bayesian Graphical Causal Models

Fit Bayesian graphical causal models specified as directed acyclic
graphs (DAGs). Each node's conditional distribution `P(X | parents)` is
modelled as a Gaussian process of its parents, exploiting the Markov
blanket factorisation. Causal effects are estimated by propagating
interventions forward through the structural model using the
do-operator, and causal structure can be discovered as a posterior
distribution over DAGs.

### Main entry point

    model <- KaguModel$new(dag = list(x = c(), y = c("x")))
    model$fit(data)
    effect <- model$effects("x", "y")
    effect$summary()

## See also

Useful links:

- <https://github.com/causalabs/kagu-r>

- <https://causalabs.github.io/kagu-r/>

- Report bugs at <https://github.com/causalabs/kagu-r/issues>

## Author

**Maintainer**: Jordan Hart <jordan@causa.tech>

Authors:

- Jordan Hart <jordan@causa.tech>

- Dan Franks <dan@causa.tech>

Other contributors:

- Causa Ltd \[copyright holder\]
