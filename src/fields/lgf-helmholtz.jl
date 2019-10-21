#using Compat.Serialization
#using Compat: @info

const GLH_NODES, GLH_WEIGHTS = gausslegendre(100)
const LGFH_DIR  = joinpath(pwd(), "cache")

alpha_to_string(α::Float64) = string(100000+α*10000)[2:6]

lgfh_file(α) = joinpath(LGFH_DIR,"lgfhtable_alpha"*alpha_to_string(α)*".dat")

function load_lgf_helmholtz(N,α)
    if isfile(lgfh_file(α))
        G = deserialize(open(lgfh_file(α), "r"))
        if size(G,1) ≥ N
            return G
        end
    end
    build_lgf_helmholtz(N,α)
end

function build_lgf_helmholtz(N,α)
    @info "Building and caching Helmholtz LGF table for α = "*string(α)

    g = zeros(N, N)
    for y in 0:N-1, x in 0:y
        g[x+1,y+1] = lgf_helmholtz(x,y)
    end

    G = Symmetric(g)
    mkpath(LGFH_DIR)
    serialize(open(lgfh_file(α),"w"),G)
    G
end

quadgauss_helmholtz(f::Function) = dot(GLH_WEIGHTS, f(GLH_NODES))

function lgf_helmholtz(i :: Integer, j :: Integer, α::Real)
  if i ≥ j
      v = ComplexF64(0.0)
      for (idex,y) in enumerate(GLH_NODES)
          xi1 = 0.5*π*(y+1)
          expixi1 = cos(i*xi1)
          cosxi1 = cos(xi1)
          v += GLH_WEIGHTS[idex]*quadgauss_helmholtz() do x
              xi2 = 0.5*π*(x .+ 1)
              expixi2 = cos.(j*xi2)
              return expixi2.*expixi1./(im*α .+ 4 .- 2*cos.(xi2) .- 2*cosxi1)
          end
      end
      return v/4
  else
      return lgf_helmholtz(j,i,α)
  end

end
