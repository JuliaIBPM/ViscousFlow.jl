using Documenter, ViscousFlow

makedocs(
    sitename = "ViscousFlow.jl",
    doctest = true,
    clean = true,
    pages = [
        "Home" => "index.md",
        "Manual" => ["manual/fields.md",
                     "manual/bodies.md",
                     "manual/saddlesystems.md",
                     "manual/timemarching.md",
                     "manual/navierstokes.md"
                     ]
        #"Internals" => [ "internals/properties.md"]
    ],
    Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
    #assets = ["assets/custom.css"],
    #strict = true
)


#if "DOCUMENTER_KEY" in keys(ENV)
deploydocs(
     repo = "github.com/jdeldre/ViscousFlow.jl.git",
     target = "build",
     deps = nothing,
     make = nothing
     #versions = "v^"
)
#end
