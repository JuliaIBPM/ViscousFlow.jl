using JLD
import JLD: save, load

export WritePlan,StorePlan, initialize_storage, store_data!


"""
    WritePlan(file,write_Δt,varlist)

Create a plan for writing data to file. The filename `file` is specified,
to be written to every `write_Δt` time units. The variable names to be written are
specified as a vector of strings in `varlist`.
"""
struct WritePlan
    filen::String
    write_Δt::Float64
    varlist::Vector{String}
end

"""
    save(t,R::WritePlan,v)

Check if time `t` is appropriate for writing to file, according to the WritePlan `R`,
and if so, write the specified variables `v` (which may be separated by commas)
in `R` to the file in `R`.
"""
function save(t,R::WritePlan,v...)
    tol = 1e-8
    if (isapprox(mod(t,R.write_Δt),0,atol=tol) ||
        isapprox(mod(t,R.write_Δt),R.write_Δt,atol=tol))
        save(R.filen,_writelist(R,v)...)
    end
end

"""
    load(R::WritePlan)

Load the data stored in the file specified by the WritePlan `R`.
"""
load(R::WritePlan) = load(R.filen)

function _writelist(R::WritePlan,v)
    list = ()
    @assert length(R.varlist)==length(v)
    for (i,vname) in enumerate(R.varlist)
        list = (list...,vname,v[i])
    end
    return list
end


"""
    StorePlan(min_t,max_t,store_Δt,v)

Create a plan for storing history data. The storage of data is specified to
start at time unit `min_t` and to proceed until (and including) `max_t`, and
is stored every `store_Δt` time units. The list of variables to be stored is
specified as a list of variables `v`. Tuple-type variables are unwrapped
into separate storage.
"""
struct StorePlan
    min_t::Float64
    max_t::Float64
    store_Δt::Float64
    varlist::Vector{DataType}
    StorePlan(min_t,max_t,store_Δt,varlist...) = new(min_t,max_t,store_Δt,_get_type(varlist))
end

#function StorePlan(min_t,max_t,store_Δt,v...)
#    return StorePlan(min_t,max_t,store_Δt,_get_type(v))
#end

"""
    initialize_storage(S::StorePlan) -> Vector

Initialize a storage data stack for the storage plan `S`. The output is
an empty vector of vectors.
"""
function initialize_storage(S::StorePlan)
    data = []
    for v in S.varlist
        push!(data,v[])
    end
    return data
end


"""
    store_data!(data,t,S::StorePlan,v)

Check whether time `t` is a time for saving for storage as described by plan `S`,
and if so, push the variables specified in `v` onto the `data` stack.
"""
function store_data!(data,t,S::StorePlan,v...)
  tol = 1e-8
  if t >= (S.min_t-tol) && t <= (S.max_t + tol) &&
      ((isapprox(mod(t,S.store_Δt),0,atol=tol) ||
        isapprox(mod(t,S.store_Δt),S.store_Δt,atol=tol)))
        store_data!(data,S::StorePlan,v)
  end
  return data
end


function _get_type(u)
    tlist = DataType[]
    _get_type!(tlist,u)
    return tlist
end
function _get_type!(tlist,u::Tuple)
    for ui in u
        _get_type!(tlist,ui)
    end
    return tlist
end
function _get_type!(tlist,u)
    push!(tlist,typeof(u))
    return tlist
end

function store_data!(data,S::StorePlan,u)
    cnt = [0]
    return _store_data!(data,cnt,S,u)
end
function _store_data!(data,cnt,S::StorePlan,u::Tuple)
    for ui in u
        _store_data!(data,cnt,S,ui)
    end
    return data
end
function _store_data!(data,cnt,S::StorePlan,u)
    cnt[1] += 1
    @assert S.varlist[cnt[1]] == typeof(u)
    push!(data[cnt[1]],u)
    return data
end
