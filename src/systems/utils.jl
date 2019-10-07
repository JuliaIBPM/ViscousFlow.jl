using JLD
import JLD: save, load

export WritePlan,StorePlan, initialize_storage, store_data!

"""
    WritePlan(file,write_Δt,varlist)

Create a plan for writing data to file. The filename `file` is specified,
to be written to every `write_Δt` time units. The variables to be written are
specified as a vector of strings in `varlist`.
"""
struct WritePlan
    filen::String
    write_Δt::Float64
    varlist::Vector{String}
end

"""
    save(t,R::WritePlan)

Check if time `t` is appropriate for writing to file, according to the WritePlan `R`,
and if so, write the specified variables in `R` to the file in `R`.
"""
function save(t,R::WritePlan)
    if (isapprox(mod(t,R.restart_Δt),0,atol=1e-6) ||
        isapprox(mod(t,R.restart_Δt),R.restart_Δt,atol=1e-6))
        save(R.filen,writelist(R)...)
    end
end

"""
    load(R::WritePlan)

Load the data stored in the file specified by the WritePlan `R`.
"""
load(R::WritePlan) = load(R.filen)

function _writelist(R::WritePlan)
    list = ()
    for v in R.varlist
        list = (list...,v,eval(Symbol(v)))
    end
    return list
end


"""
    StorePlan(min_t,max_t,store_Δt,varlist)

Create a plan for storing history data. The storage of data is specified to
start at time unit `min_t` and to proceed until (and including) `max_t`, and
is stored every `store_Δt` time units. The list of variables to be stored is
specified as a vector of strings in `varlist`. Tuple-type variables are unwrapped
into separate storage.
"""
struct StorePlan
    min_t::Float64
    max_t::Float64
    store_Δt::Float64
    varlist::Vector{String}
end

"""
    initialize_storage(S::StorePlan) -> Vector

Initialize a storage data stack for the storage plan `S`. The output is
an empty vector of vectors.
"""
function initialize_storage(R::StorePlan)
    data = []
    for (i,v) in enumerate(R.varlist)
        vtype = typeof(eval(Symbol(v)))
        if vtype <: Tuple
            for vi in eval(Symbol(v))
                push!(data,typeof(vi)[])
            end
        else
            push!(data,vtype[])
        end
    end
    return data
end

"""
    store_data!(data,t,S::StorePlan)

Check whether time `t` is a time for saving for storage as described by plan `S`,
and if so, push the variables specified in `S` onto the `data` stack.
"""
function store_data!(data,t,R::StorePlan)
    if t >= (R.min_t-1e-6) && t <= (R.max_t + 1e-6) &&
        ((isapprox(mod(t,R.store_Δt),0,atol=1e-6) ||
          isapprox(mod(t,R.store_Δt),R.store_Δt,atol=1e-6)))
        cnt = 0
        for v in R.varlist
            if typeof(eval(Symbol(v))) <: Tuple
                for vi in eval(Symbol(v))
                    cnt += 1
                    push!(data[cnt],vi)
                end
            else
                cnt += 1
                push!(data[cnt],eval(Symbol(v)))
            end
        end
    end
    return data
end
