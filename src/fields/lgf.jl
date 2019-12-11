using FastGaussQuadrature
using Serialization
#using Compat: @info

const GL_NODES, GL_WEIGHTS = gausslegendre(100)
const LGF_DIR  = joinpath(@__DIR__, "cache")
const LGF_FILE = joinpath(LGF_DIR, "lgftable.dat")

function load_lgf(N)
    if isfile(LGF_FILE)
        G = deserialize(open(LGF_FILE, "r"))
        if size(G,1) ≥ N
            return G
        end
    end
    build_lgf(N)
end

function build_lgf(N)
    @info "Building and caching LGF table"

    g = zeros(N, N)
    for y in 0:N-1, x in 0:y
        g[x+1,y+1] = lgf(x,y)
    end

    G = Symmetric(g)
    mkpath(LGF_DIR)
    #open(LGF_FILE, "w") do f
    #    serialize(f, G)
    #end
    serialize(open(LGF_FILE,"w"),G)
    G
end

quadgauss(f::Function) = dot(GL_WEIGHTS, f(GL_NODES))

function lgf(i, j)
    if i == j ==0
        return 0.0
    elseif i ≥ j
        v = quadgauss() do x
            if x == -1
                return sqrt(2)abs(i)
            else
                t = (x .+ 1)./2
                return 0.5real((1 .-
                                ( (t.-sqrt(1im))./(t.+sqrt(1im)) ).^(j.+abs(i)).*
                                ( (t.+sqrt(-1im))./(t.-sqrt(-1im)) ).^(j.-abs(i)) ))./t
            end

        end
        return 0.5v/pi
    else
        return lgf(j,i)
    end

end

const LGF_TABLE = load_lgf(1600)
