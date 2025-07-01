#= Fields =#

vorticity!(masked_w::Nodes{Dual},w::Nodes{Dual},dv::VectorData,base_cache::BasicILMCache,wcache::VectorPotentialCache) =
            masked_curlv_from_curlv_masked!(masked_w,w,dv,base_cache,wcache)

total_vorticity!(w::Nodes{Dual},masked_w::Nodes{Dual},dv::VectorData,base_cache::BasicILMCache,wcache::VectorPotentialCache) =
            curlv_masked_from_masked_curlv!(w,masked_w,dv,base_cache,wcache)

for f in [:vorticity,:total_vorticity]

    f! = Symbol(string(f)*"!")

    @eval function $f!(masked_w::Nodes{Dual},w::Nodes{Dual},τ,x,sys::ILMSystem,t)
        @unpack extra_cache, base_cache = sys
        @unpack dvb, velcache = extra_cache
        @unpack wcache = velcache

        prescribed_surface_jump!(dvb,x,t,sys)
        $f!(masked_w,w,dvb,base_cache,wcache)
    end

    @eval $f(w::Nodes{Dual},τ,x,sys::ILMSystem,t) = $f!(zeros_gridcurl(sys),w,τ,x,sys,t)

    @eval @snapshotoutput $f
end

function velocity!(v::Edges{Primal},curl_vmasked::Nodes{Dual},masked_divv::Nodes{Primal},dv::VectorData,vp,base_cache::BasicILMCache,velcache::VectorFieldCache,w_tmp::Nodes{Dual})
    @unpack wcache  = velcache

    masked_curlv = w_tmp

    masked_curlv_from_curlv_masked!(masked_curlv,curl_vmasked,dv,base_cache,wcache)
    vecfield_helmholtz!(v,masked_curlv,masked_divv,dv,vp,base_cache,velcache)

    return v
end

function velocity!(v::Edges{Primal},w::Nodes{Dual},x,sys::ILMSystem,t)
    @unpack phys_params, extra_cache, base_cache = sys
    @unpack dvb, velcache, divv_tmp, w_tmp = extra_cache

    prescribed_surface_jump!(dvb,x,t,sys)

    Vinf = evaluate_freestream(t,x,sys)

    fill!(divv_tmp,0.0)
    fill!(w_tmp,0.0)
    velocity!(v,w,divv_tmp,dvb,Vinf,base_cache,velcache,w_tmp)
end

velocity(w::Nodes{Dual},τ,x,sys::ILMSystem,t) = velocity!(zeros_grid(sys),w,x,sys,t)

@snapshotoutput velocity


function streamfunction!(ψ::Nodes{Dual},w::Nodes{Dual},vp::Tuple,base_cache::BasicILMCache,wcache::VectorPotentialCache)
    @unpack stemp = wcache

    fill!(ψ,0.0)
    vectorpotential_from_curlv!(stemp,w,base_cache)
    ψ .+= stemp

    vectorpotential_uniformvecfield!(stemp,vp...,base_cache)
    ψ .+= stemp

    return ψ
end

function streamfunction!(ψ::Nodes{Dual},w::Nodes{Dual},x,sys::ILMSystem,t)
    @unpack phys_params, extra_cache, base_cache = sys
    @unpack velcache = extra_cache
    @unpack wcache = velcache

    Vinf = evaluate_freestream(t,x,sys)

    streamfunction!(ψ,w,Vinf,base_cache,wcache)

end

streamfunction(w::Nodes{Dual},τ,x,sys::ILMSystem,t) = streamfunction!(zeros_gridcurl(sys),w,x,sys,t)

@snapshotoutput streamfunction

function pressure!(press::Nodes{Primal},w::Nodes{Dual},τ,x,sys::ILMSystem,t)
      @unpack extra_cache, base_cache, phys_params, motions = sys
      @unpack velcache, v_tmp, dv_tmp, dv, divv_tmp, v_rot = extra_cache
      @unpack reference_body, vl, Xl, m = motions

      velocity!(v_tmp,w,x,sys,t)
      viscousflow_velocity_ode_rhs!(dv,v_tmp,x,sys,t)
      fill!(dv_tmp,0.0)
      viscousflow_velocity_constraint_force!(dv_tmp,τ,x,sys)
      dv .-= dv_tmp
      viscousflow_velocity_ode_implicit_rhs!(dv_tmp,x,sys,t)
      dv .-= dv_tmp

      divergence!(divv_tmp,dv,base_cache)
      inverse_laplacian!(press,divv_tmp,base_cache)

      # compute R^T v = v' - R^T(Ur - Uinf) here (velocity field in inertial frame,
      #    but in rotating coordinates)
      #velocity_rel_to_inertial_frame!(v_tmp,x,t,base_cache,phys_params,motions)
      press .-= 0.5*magsq(v_tmp)

      # For rotating coordinate systems...
      if reference_body != 0
        velocity_rel_to_inertial_frame!(v_tmp,x,t,base_cache,phys_params,motions)        

        # Compute Ω × x̂ here
        fill!(v_rot,0.0)
        velocity_rel_to_rotating_frame!(v_rot,x,t,base_cache,phys_params,motions)
        v_rot .*= -1.0

        # Now this holds rigid body field, V̂r = Ω × x̂ + R^T(Ur - Uinf)
        velocity_rel_to_inertial_frame!(v_rot,x,t,base_cache,phys_params,motions)

        # Now compute V̂r ⋅ R^T v
        fill!(divv_tmp,0.0)
        gridwise_dot!(divv_tmp,v_rot,v_tmp)

        #press .-= 0.5*magsq(v_tmp) - divv_tmp
        press .+= divv_tmp



        #=
        # Compute v̂ here
        velocity_rel_to_rotating_frame!(v_tmp,x,t,base_cache,phys_params,motions)

        press .-= 0.5*magsq(v_tmp) - 0.5*magsq(v_rot)

        # Compute U̇r⋅x̂ here. This is actually (Ω x x̂)⋅U_r
        vlref = motions.vl[reference_body]
        vx, vy = vlref.linear
        v_rot.u .*= vx
        v_rot.v .*= vy
        fill!(divv_tmp,0.0)
        grid_interpolate!(divv_tmp,v_rot)
        press .-= divv_tmp
        =#

      end

      return press
end

pressure(w::Nodes{Dual},τ,x,sys::ILMSystem,t) = pressure!(zeros_griddiv(sys),w,τ,x,sys,t)

@snapshotoutput pressure

function convective_acceleration!(vdv::Edges{Primal},w::Nodes{Dual},x,sys::ILMSystem,t)
    @unpack extra_cache, base_cache = sys
    @unpack v_tmp, cdcache = extra_cache
    velocity!(v_tmp,w,x,sys,t)
    vdv .= convective_derivative(v_tmp,base_cache)
    return vdv
end

convective_acceleration(w::Nodes{Dual},τ,x,sys::ILMSystem,t) = convective_acceleration!(zeros_grid(sys),w,x,sys,t)

@snapshotoutput convective_acceleration


function Qcrit!(Q::Nodes{Primal},w::Nodes{Dual},τ,x,sys::ILMSystem,t)
    @unpack motions = sys
    vel = zeros_grid(sys)

    velocity!(vel,w,x,sys,t)

    gradv = zeros_gridgrad(sys)
    grad!(gradv,vel,sys)
    tmp = zeros_griddiv(sys)
    grid_interpolate!(tmp,2*gradv.dudy∘gradv.dvdx)

    Q .= gradv.dudx∘gradv.dudx + gradv.dvdy∘gradv.dvdy + tmp
    Q .*= -0.5
    return Q
end

Qcrit(w::Nodes{Dual},τ,x,sys::ILMSystem,t) = Qcrit!(zeros_griddiv(sys),w,τ,x,sys,t)

@snapshotoutput Qcrit

#= Surface fields =#

function traction!(tract::VectorData{N},τ::VectorData{N},x,sys::ILMSystem,t) where {N}
    @unpack bc, extra_cache, base_cache, phys_params, motions = sys
    @unpack vb_tmp, dvb = extra_cache
    @unpack sscalar_cache, bl = base_cache
    @unpack reference_body, m = motions

    # (v̅ - Ẋ)⋅n -> sscalar_cache
    prescribed_surface_average!(vb_tmp,x,t,sys)
    if reference_body != 0
        surface_velocity_in_translating_frame!(dvb,x,base_cache,motions,t)
        vb_tmp .-= dvb
        velocity_rel_to_rotating_frame!(vb_tmp,x,t,base_cache,phys_params,motions)
    else
        surface_velocity!(dvb,x,sys,t)
        vb_tmp .-= dvb
    end

    bl_tmp = deepcopy(bl)
    update_body!(bl_tmp,x,m)

    nrm = normals(bl_tmp)
    pointwise_dot!(sscalar_cache,nrm,vb_tmp)

    # [v] -> dvb
    prescribed_surface_jump!(dvb,x,t,sys)

    # [v](v̅ - Ẋ)⋅n
    product!(tract,dvb,sscalar_cache)

    tract .+= τ

    return tract

end
traction!(out::VectorData{0},τ::VectorData{0},x,sys::ILMSystem,t) = out
traction(w::Nodes{Dual},τ::VectorData,x,sys::ILMSystem,t) = traction!(zeros_surface(sys),τ,x,sys,t)
@snapshotoutput traction

function pressurejump!(dpb::ScalarData{N},τ::VectorData{N},x,sys::ILMSystem,t) where {N}
    @unpack base_cache, motions = sys
    @unpack sdata_cache, bl = base_cache
    @unpack m = motions

    fill!(dpb,0.0)
    fill!(sdata_cache,0.0)
    bl_tmp = deepcopy(bl)
    update_body!(bl_tmp,x,m)

    nrm = normals(bl_tmp)
    traction!(sdata_cache,τ,x,sys,t)
    pointwise_dot!(dpb,nrm,sdata_cache)
    dpb .*= -1.0
    return dpb
end
pressurejump!(out::VectorData{0},τ::VectorData{0},x,sys::ILMSystem,t) = out
pressurejump(w::Nodes{Dual},τ::VectorData,x,sys::ILMSystem,t) = pressurejump!(zeros_surfacescalar(sys),τ,x,sys,t)
@snapshotoutput pressurejump

function pressureplus(w::Nodes{Dual},τ::VectorData,x,sys::ILMSystem,t)
    @unpack base_cache = sys
    @unpack Ediv = base_cache
    pbarb = Ediv*pressure(w,τ,x,sys,t)
    dpb = pressurejump(w,τ,x,sys,t)
    return pbarb + 0.5*dpb
end
@snapshotoutput pressureplus

function pressureminus(w::Nodes{Dual},τ::VectorData,x,sys::ILMSystem,t)
    @unpack base_cache = sys
    @unpack Ediv = base_cache
    pbarb = Ediv*pressure(w,τ,x,sys,t)
    dpb = pressurejump(w,τ,x,sys,t)
    return pbarb - 0.5*dpb
end
@snapshotoutput pressureminus


function shearstressjump!(dtaub::ScalarData{N},τ::VectorData{N},x,sys::ILMSystem,t) where {N}
    @unpack base_cache, motions = sys
    @unpack sdata_cache, bl = base_cache
    @unpack m = motions

    bl_tmp = deepcopy(bl)
    update_body!(bl_tmp,x,m)

    nrm = normals(bl_tmp)
    traction!(sdata_cache,τ,x,sys,t)
    pointwise_cross!(dtaub,nrm,sdata_cache)
    return dtaub
end
shearstressjump!(out::VectorData{0},τ::VectorData{0},x,sys::ILMSystem,t) = out
shearstressjump(w::Nodes{Dual},τ::VectorData,x,sys::ILMSystem,t) = shearstressjump!(zeros_surfacescalar(sys),τ,x,sys,t)
@snapshotoutput shearstressjump

#= Integrated metrics =#

force(w::Nodes{Dual},τ::VectorData{0},x,sys::ILMSystem{S,P,0},t,bodyi::Int;kwargs...) where {S,P} = nothing, nothing #Vector{Float64}(), Vector{Float64}()

"""
    force(sol,sys,bodyi[;axes=0,reference_body=bodyi]) -> Tuple{Vector}

Calculate the moment and force exerted by the fluid on body `bodyi` from the computational solution `sol` of system `sys`.
It returns the force history as a tuple of arrays: one array for each component.
If `axes=0` (the default), then the components are provided in the inertial coordinate system. However,
this can be changed to any body ID to provide the components in another system. The
keyword `force_reference` determines which body's system the moment is taken about.
By default, it is the body `bodyi` itself, but any other body can be specified,
or the inertial system (by specifying 0).
""" force(sol,sys,bodyi)

function force(w::Nodes{Dual},τ::VectorData{N},x,sys::ILMSystem{S,P,N},t,bodyi::Int;axes=0,force_reference=bodyi,force_type=:all) where {S,P,N}
    @unpack base_cache, extra_cache, phys_params, motions = sys
    @unpack sdata_cache, sscalar_cache, bl = base_cache
    @unpack vb_tmp = extra_cache
    @unpack reference_body, m, Xl = motions

    bl_tmp = deepcopy(bl)
    update_body!(bl_tmp,x,m)
    ds = areas(bl_tmp)

    _force_integrand!(sdata_cache,w,τ,x,sys,t,sscalar_cache,Val(force_type))
    fx = integrate(sdata_cache.u,ds,bl_tmp,bodyi)
    fy = integrate(sdata_cache.v,ds,bl_tmp,bodyi)

    vb_tmp .= deepcopy(points(bl_tmp))
    xc, yc = bl_tmp[bodyi].cent
    vb_tmp.u .-= xc
    vb_tmp.v .-= yc

    pointwise_cross!(sscalar_cache,vb_tmp,sdata_cache)
    mom = integrate(sscalar_cache,ds,bl_tmp,bodyi)

    if !(force_reference == bodyi && axes == reference_body)

        fb = PluckerForce{2}(angular=mom,linear=[fx,fy])
        evaluate_motion!(motions,x,t)
        Xfb_to_bref = force_transform_from_A_to_B(Xl,bodyi,force_reference)

        Rfbref_to_axes = rotation_transform(force_transform_from_A_to_B(Xl,force_reference,axes))

        fbref = Rfbref_to_axes*Xfb_to_bref*fb
        mom, fx, fy = fbref
    end

    return mom, fx, fy
end

_force_integrand!(data,w,τ,x,sys,t,scalardata,::Val{:all}) = traction!(data,τ,x,sys,t)

function _force_integrand!(data,w,τ,x,sys,t,scalardata,::Val{:pressure})
    nrm = normals(sys)
    #scalardata .= pressureplus(w,τ,x,sys,t)
    pressurejump!(scalardata,τ,x,sys,t)
    product!(data,-nrm,scalardata)
    return nothing
end

function _force_integrand!(data,w,τ,x,sys,t,scalardata,::Val{:viscous})
    nrm = normals(sys)
    shearstressjump!(scalardata,τ,x,sys,t)
    pointwise_cross!(data,-nrm,scalardata)
    return nothing
end

@vectorsurfacemetric force

moment(w::Nodes{Dual},τ::VectorData{0},x,sys::ILMSystem{S,P,0},t,bodyi::Int;kwargs...) where {S,P} = nothing # Vector{Float64}()

"""
    moment(sol,sys,bodyi[;center=(0,0)]) -> Vector

Calculated the momemt exerted by the fluid on body `bodyi` from the computational solution `sol` of system `sys`.
It returns the moment history as an array. The moment is calculated about center
`center`, which defaults to `(0,0)` in whichever coordinate system the problem is solved.
""" moment(sol,sys,bodyi)


function moment(w::Nodes{Dual},τ::VectorData{N},x,sys::ILMSystem{S,P,N},t,bodyi::Int;center=(0.0,0.0)) where {S,P,N}
    @unpack base_cache, extra_cache = sys
    @unpack sdata_cache, sscalar_cache, bl = base_cache
    @unpack vb_tmp = extra_cache

    bl_tmp = deepcopy(bl)
    update_body!(bl_tmp,x,m)
    ds = areas(bl_tmp)

    xc, yc = center
    vb_tmp .= deepcopy(points(sys))
    vb_tmp.u .-= xc
    vb_tmp.v .-= yc

    traction!(sdata_cache,τ,x,sys,t)
    pointwise_cross!(sscalar_cache,vb_tmp,sdata_cache)
    mom = integrate(sscalar_cache,ds,bl_tmp,bodyi)

    return mom
end

@scalarsurfacemetric moment


power(w::Nodes{Dual},τ::VectorData{0},x,sys::ILMSystem{S,P,0},t,bodyi::Int;kwargs...) where {S,P} = nothing # Vector{Float64}()

"""
    power(sol,sys,bodyi;include_freestream=false)

Calculate the history of the total rate of work done by the flow on body `bodyi` (or on the flow by the body, if negative)
from the computational solution `sol` of system `sys`. If `freestream=false` (default), then
the free stream is not subtracted from the translational velocity in the calculation. If
it is `true`, then the translational part of power is based on translational velocity
minus the free stream (i.e., a reference frame at rest at infinity).
""" power(sol,sys,bodyi)


function power(w::Nodes{Dual},τ::VectorData{N},x,sys::ILMSystem{S,P,N},t,bodyi::Int; include_freestream=false) where {S,P,N}
    @unpack phys_params, motions = sys
    @unpack m = motions

    # Get the force and moment in the body's own coordinate system
    f = PluckerForce([force(w,τ,x,sys,t,bodyi;axes=bodyi,force_reference=bodyi)...])

    # Get velocities in the body's own coordinate system
    evaluate_motion!(motions,x,t)
    v = motions.vl[bodyi]
  
    # subtract freestream from the translational velocity
    if include_freestream
        freestream_func = get_freestream_func(phys_params)
        XRref = rotation_transform(motions.Xl[bodyi])
        Vinf_plucker = XRref*PluckerMotion{2}(linear = [freestream_func(t,phys_params)...])
        pow = dot(f,v-Vinf_plucker)
    else
        pow = dot(f,v)
    end

    return pow

end

@scalarsurfacemetric power

