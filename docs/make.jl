using Documenter, ViscousFlow

ENV["GKSwstype"] = "nul" # removes GKS warnings during plotting

makedocs(
    sitename = "ViscousFlow.jl",
    doctest = true,
    clean = true,
    pages = [
        "Home" => "index.md",
        "Manual" => ["manual/1.-Basic-viscous-flow.md",
                     "manual/2.-Basic-flow-with-a-stationary-body.md",
                     "manual/3.-Applying-pulse-forcing-to-a-flow.md",
                     "manual/4.-Multiple-stationary-bodies.md",
                     "manual/5.-Viscous-flow-about-a-moving-body.md",
                     "manual/6.-Variable-free-stream.md",
                     "manual/functions.md"
                     ]
        #"Internals" => [ "internals/properties.md"]
    ],
    #format = Documenter.HTML(assets = ["assets/custom.css"])
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict()
            )
        ))
    ),
    #assets = ["assets/custom.css"],
    #strict = true
)


#if "DOCUMENTER_KEY" in keys(ENV)
deploydocs(
     repo = "github.com/JuliaIBPM/ViscousFlow.jl.git",
     target = "build",
     deps = nothing,
     make = nothing
     #versions = "v^"
)
#end
