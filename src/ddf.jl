module DDF

export ddf_roma

roma1(r) = (1+sqrt(-3r^2+1))/3;
roma2(r) = (5-3r-sqrt(1-3(1-r)^2))/6;

function ddf_roma(r::Float64)
    rr = abs(r)
    if rr <= 0.5
        roma1(rr)
    elseif 0.5 < rr < 1.5
    	roma2(rr)
    else
        0.0
    end
end

function ddf_roma(r::Array{Float64,1})
  rr = abs.(r)

  f = zeros(rr)
  pts = find(x -> x<=0.5,rr)
  f[pts] = roma1.(rr[pts])

  pts = find(x -> 0.5 < x < 1.5,rr)
  f[pts] = roma2.(rr[pts])

  f

end

end
