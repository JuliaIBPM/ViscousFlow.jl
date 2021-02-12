### Utilities ###

struct PulseParams{GF}
    field :: GF
    t0 :: Float64
    σt :: Float64
end

_process_pulses(::Nothing,s,grid) = nothing

function _process_pulses(params::PulseParams,s,grid)
  gf = GeneratedField(s,params.field,grid)
  pulse = PulseField(gf,params.t0,params.σt)
  return [pulse]
end

function _process_pulses(params::Vector{<:PulseParams},s,grid)
  pulses = ModulatedField[]
  for p in params
    gf = GeneratedField(s,p.field,grid)
    push!(pulses,PulseField(gf,p.t0,p.σt))
  end
  return pulses
end
