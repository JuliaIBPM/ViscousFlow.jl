using ViscousFlow
using Test
##using TestSetExtensions
using Literate
using Suppressor

const GROUP = get(ENV, "GROUP", "All")

notebookdir = "../examples"
docdir = "../docs/src/manual"
litdir = "./literate"

if GROUP == "All" || GROUP == "Auxiliary"
  include("pointforce.jl")
end


if GROUP == "All" || GROUP == "Notebooks"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      endswith(file,".jl") && startswith(file,"3") && Literate.notebook(joinpath(root, file),notebookdir)
      #endswith(file,".jl") && Literate.notebook(joinpath(root, file),notebookdir)
    end
  end
end

if GROUP == "Documentation"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      endswith(file,".jl") && Literate.markdown(joinpath(root, file),docdir)
    end
  end
end
