# Enumerate all DAGs over a set of nodes

Generates every directed acyclic graph over `nodes`, optionally
excluding a set of disallowed directed edges. Each returned DAG is in
the named-list format used throughout the package (node -\> character
vector of parents).

## Usage

``` r
enumerate_dags(nodes, disallowed = NULL, required = NULL, allow_empty = FALSE)
```

## Arguments

- nodes:

  Character vector of node names.

- disallowed:

  Optional list of length-2 character vectors `c(from, to)`, each
  forbidding the directed edge `from -> to`. Directions are independent:
  forbidding `a -> b` still permits `b -> a`.

- required:

  Optional list of length-2 character vectors `c(from, to)`, each
  requiring the directed edge `from -> to` to be present in all
  candidate DAGs.

- allow_empty:

  Logical - whether to include the completely empty (edgeless) graph in
  the search space. Defaults to `FALSE`.

## Value

A list of DAG specifications. Over 3 unconstrained nodes there are 24
(with `allow_empty = FALSE`).

## Details

The number of DAGs grows super-exponentially in the number of nodes, so
the search space is only tractable for small problems. Use `disallowed`
(e.g. to encode a known temporal ordering) to prune it.

## Examples

``` r
length(enumerate_dags(c("a", "b", "c")))            # 24
#> [1] 24
# forbid b -> a and c -> a (e.g. 'a' comes first in time)
enumerate_dags(c("a", "b", "c"),
               disallowed = list(c("b", "a"), c("c", "a")))
#> [[1]]
#> [[1]]$a
#> character(0)
#> 
#> [[1]]$b
#> [1] "a"
#> 
#> [[1]]$c
#> character(0)
#> 
#> 
#> [[2]]
#> [[2]]$a
#> character(0)
#> 
#> [[2]]$b
#> [1] "c"
#> 
#> [[2]]$c
#> character(0)
#> 
#> 
#> [[3]]
#> [[3]]$a
#> character(0)
#> 
#> [[3]]$b
#> [1] "a" "c"
#> 
#> [[3]]$c
#> character(0)
#> 
#> 
#> [[4]]
#> [[4]]$a
#> character(0)
#> 
#> [[4]]$b
#> character(0)
#> 
#> [[4]]$c
#> [1] "a"
#> 
#> 
#> [[5]]
#> [[5]]$a
#> character(0)
#> 
#> [[5]]$b
#> [1] "a"
#> 
#> [[5]]$c
#> [1] "a"
#> 
#> 
#> [[6]]
#> [[6]]$a
#> character(0)
#> 
#> [[6]]$b
#> [1] "c"
#> 
#> [[6]]$c
#> [1] "a"
#> 
#> 
#> [[7]]
#> [[7]]$a
#> character(0)
#> 
#> [[7]]$b
#> [1] "a" "c"
#> 
#> [[7]]$c
#> [1] "a"
#> 
#> 
#> [[8]]
#> [[8]]$a
#> character(0)
#> 
#> [[8]]$b
#> character(0)
#> 
#> [[8]]$c
#> [1] "b"
#> 
#> 
#> [[9]]
#> [[9]]$a
#> character(0)
#> 
#> [[9]]$b
#> [1] "a"
#> 
#> [[9]]$c
#> [1] "b"
#> 
#> 
#> [[10]]
#> [[10]]$a
#> character(0)
#> 
#> [[10]]$b
#> character(0)
#> 
#> [[10]]$c
#> [1] "a" "b"
#> 
#> 
#> [[11]]
#> [[11]]$a
#> character(0)
#> 
#> [[11]]$b
#> [1] "a"
#> 
#> [[11]]$c
#> [1] "a" "b"
#> 
#> 
```
