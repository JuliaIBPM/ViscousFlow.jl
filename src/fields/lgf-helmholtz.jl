#using Compat.Serialization
#using Compat: @info

# Sets up a variety of different Gauss-Legendre data lengths
# for resolving different regimes of oscillatory function
#const GLH_NODES_20, GLH_WEIGHTS_20 = gausslegendre(20)
const GLH_NODES_100, GLH_WEIGHTS_100 = gausslegendre(100)
const GLH_NODES_200, GLH_WEIGHTS_200 = gausslegendre(200)
const GLH_NODES_500, GLH_WEIGHTS_500 = gausslegendre(500)
const GLH_NODES_1000, GLH_WEIGHTS_1000 = gausslegendre(1000)
const GLH_NODES_2000, GLH_WEIGHTS_2000 = gausslegendre(2000)
const GLH_NODES_2500, GLH_WEIGHTS_2500 = gausslegendre(2500)
const GLH_NODES_3500, GLH_WEIGHTS_3500 = gausslegendre(3500)
const GLH_NODES_10000, GLH_WEIGHTS_10000 = gausslegendre(10000)
const GLH_N = [100,200,500,1000,2000,2500,3500,10000]

const LGFH_DIR  = joinpath(pwd(), "cache")

using ProgressMeter

alpha_to_string(α::Float64) = string(100000+α*10000)[2:6]

lgfh_file(α) = joinpath(LGFH_DIR,"lgfhtable_alpha"*alpha_to_string(α)*".dat")

function _get_gl_data(n;i=0)
  # Set the number of Gauss points so that it will resolve the Nyquist frequency
  # on the oscillatory part of the Fourier integral
  while 2*n > GLH_N[i += 1] end
  return getfield(Fields,Symbol("GLH_NODES_",string(GLH_N[i]))), getfield(Fields,Symbol("GLH_WEIGHTS_",string(GLH_N[i])))
end

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


    asymptotic_radius = _get_lgf_switch(α)
    @info "Switch to asymptotic formula at i = "*string(asymptotic_radius)

    lgf_support = _get_lgf_support(α)
    @info "Support radius = "*string(lgf_support)


    g = zeros(ComplexF64,N, N)
    @showprogress for y in 0:N-1, x in 0:y
      ir = _index_radius(x,y)
      if ir < asymptotic_radius
        g[x+1,y+1] = lgf_helmholtz(x,y,α)
      elseif ir < lgf_support
        g[x+1,y+1] = lgf_helmholtz_asymptotic(x,y,α)
      end
    end

    G = Symmetric(g)
    mkpath(LGFH_DIR)
    serialize(open(lgfh_file(α),"w"),G)
    G
end

#quadgauss_helmholtz(f::Function) = dot(GLH_WEIGHTS, f(GLH_NODES))

function lgf_helmholtz(i :: Integer, j :: Integer, α::Real)

  nodes_i, weights_i = _get_gl_data(abs(i))
  nodes_j, weights_j = _get_gl_data(abs(j))

  quadgauss_j(f::Function) = dot(weights_j, f(nodes_j))

  if i ≥ j
      v = ComplexF64(0.0)
      for (idex,y) in enumerate(nodes_i)
          xi1 = 0.5*π*(y+1)
          expixi1 = cos(i*xi1)
          cosxi1 = cos(xi1)
          v += weights_i[idex]*quadgauss_j() do x
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

lgf_helmholtz_asymptotic(i :: Integer, j :: Integer, α::Real) =
      conj(im*0.25*hankelh1(0,exp(im*π/4)*sqrt(α)*sqrt(i^2+j^2)))

function _get_lgf_support(α::Number;tol=eps(Float64))
  i = 0
  while abs(lgf_helmholtz(i,0,α)) > tol
     i += 1
  end
  return i
end

function _get_lgf_switch(α::Number;tol=1e-9)
  i = 1
  while abs(lgf_helmholtz(i,0,α)-lgf_helmholtz_asymptotic(i,0,α)) > tol
     i += 1
  end
  return i
end

_index_radius(i::Int,j::Int) = ceil(Int,sqrt(i^2+j^2))
