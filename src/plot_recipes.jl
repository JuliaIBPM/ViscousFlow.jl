using RecipesBase
using ColorTypes


@recipe function f(w::T,sys::NavierStokes) where {T<:GridData}
  @series begin
    trim := 2
    w, sys.grid
  end

  if !(isnothing(sys.bodies))
    if typeof(w) <: ScalarGridData
      @series begin
        sys.bodies
      end
    else
      if typeof(w) <: VectorGridData
        nsub = 2
      elseif typeof(w) <: TensorGridData
        nsub = 4
      end
      for i in 1:nsub
        @series begin
          subplot := i
          sys.bodies
        end
      end
    end
  end
end

@recipe function f(w1::T1,w2::T2,sys::NavierStokes) where {T1<:GridData,T2<:GridData}
    @series begin
      w1, sys.grid
    end

    @series begin
      linestyle --> :dash
      w2, sys.grid
    end

    if !(isnothing(sys.bodies))
      @series begin
        sys.bodies
      end
    end
end

@recipe f(kin::Kinematics) = RigidBodyMotion(kin)
