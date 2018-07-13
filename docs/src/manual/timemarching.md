# Time marching

```@meta
DocTestSetup = quote
using Whirl
end
```

```math
\def\ddt#1{\frac{\mathrm{d}#1}{\mathrm{d}t}}

\renewcommand{\vec}{\boldsymbol}
\newcommand{\uvec}[1]{\vec{\hat{#1}}}
\newcommand{\utangent}{\uvec{\tau}}
\newcommand{\unormal}{\uvec{n}}

\renewcommand{\d}{\,\mathrm{d}}
```


```@setup create
using Whirl
using Plots
```


## Methods

```@autodocs
Modules = [TimeMarching]
Order   = [:type, :function]
```

## Index

```@index
Pages = ["timemarching.md"]
```
