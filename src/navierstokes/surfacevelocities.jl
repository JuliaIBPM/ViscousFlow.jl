# Routines to compute surface velocities and their jumps

@inline function surface_velocity_jump!(Δus::VectorData{N},sys::NavierStokes{NX,NY,N},t::Real) where {NX,NY,N}
    @unpack surfacevel_cache, motions, bodies = sys
    @unpack Vb = surfacevel_cache
    assign_velocity!(Vb,bodies,motions,t)
    surface_velocity_jump!(Δus,Vb,sys)
    return Δus
  end

@inline surface_velocity_jump!(Δus::VectorData{N},us::VectorData{N},
                               ::NavierStokes{NX,NY,N,MT,FS,InternalFlow}) where {NX,NY,N,MT,FS} =
                               (Δus .= -us; Δus)

@inline surface_velocity_jump!(Δus::VectorData{N},us::VectorData{N},
                               ::NavierStokes{NX,NY,N,MT,FS,ExternalFlow}) where {NX,NY,N,MT,FS} =
                               (Δus .= us; Δus)

@inline surface_velocity_jump!(Δus::VectorData{N},us::VectorData{N},
                               ::NavierStokes{NX,NY,N,MT,FS,ExternalInternalFlow}) where {NX,NY,N,MT,FS} =
                               (Δus .= 0.0; Δus)


@inline function surface_velocity!(ūs::VectorData{N},sys::NavierStokes{NX,NY,N},t::Real) where {NX,NY,N}
    @unpack surfacevel_cache, motions, bodies = sys
    @unpack Vb = surfacevel_cache
    assign_velocity!(Vb,bodies,motions,t)
    surface_velocity!(ūs,Vb,sys)
    return ūs
  end

@inline surface_velocity!(ūs::VectorData{N},us::VectorData{N},
                               ::NavierStokes{NX,NY,N,MT,FS,ExternalInternalFlow}) where {NX,NY,N,MT,FS} =
                               (ūs .= us; ūs)

@inline surface_velocity!(ūs::VectorData{N},us::VectorData{N},
                               ::NavierStokes{NX,NY,N,MT,FS,SD}) where {NX,NY,N,MT,FS,SD} =
                               (ūs .= 0.5*us; ūs)

@inline function relative_surface_velocity!(ūsr::VectorData{N},sys::NavierStokes{NX,NY,N},t::Real) where {NX,NY,N}
  @unpack surfacevel_cache, motions, bodies = sys
  @unpack Vb = surfacevel_cache
  assign_velocity!(Vb,bodies,motions,t)
  surface_velocity!(ūsr,Vb,sys)
  ūsr .-= Vb
end
