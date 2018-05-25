using Documenter, Whirl

makedocs(
    format =:html,
    sitename = "Whirl.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => ["manual/fields.md"
        #             "manual/elements.md",
        #             "manual/velocities.md",
        #             "manual/timemarching.md",
        #             "manual/noflowthrough.md",
        #             "manual/motions.md"
                     ]
        #"Internals" => [ "internals/properties.md"]
    ],
    assets = ["assets/custom.css"],
    strict = true
)


if "DOCUMENTER_KEY" in keys(ENV)
    deploydocs(
     repo = "github.com/jdeldre/Whirl.jl.git",
     target = "build",
     deps = nothing,
     make = nothing,
     julia = "0.6"
    )
end
