
function tMPS(C, h::Vector{Array{Float64,2}}, t::Float64, Δt::Float64, D::Int64=Inf)
    N=div(t,Δt);
    d=size(C[1])[3];
    L=size(C)[1]
    Op_1=expm(h[2]*Δt/2.)
    Op_2=expm(h[2]*Δt)
    Op_1_1=expm(h[1]*Δt/2.)
    Op_1_l=expm(h[3]*Δt/2.)
    Op_2_l=expm(h[3]*Δt)


    asym_operator_2bond(C,:uneven,Op_1,Op_1_1,Op_1_l)
    #as Op_1 commutes with itself, if N>1 two consecutive steps with Op_1 can be combined
    for i=1:N-1
      compress(C,D,:R;stochastic=true)
      compress(C,D,:L;stochastic=true)
      asym_operator_2bond(C,:even,Op_2,Op_2_l)
      compress(C,D,:L;stochastic=true)
      compress(C,D,:R;stochastic=true)
      asym_operator_2bond(C,:uneven,Op_1,Op_1_1,Op_1_l)
    end
    compress(C,D,:L;stochastic=true)
    compress(C,D,:R;stochastic=true)
    asym_operator_2bond(C,:even,Op_2,Op_2_l)
    compress(C,D,:R;stochastic=true)
    compress(C,D,:L;stochastic=true)
    asym_operator_2bond(C,:uneven,Op_1,Op_1_1,Op_1_l)
    compress(C,D,:R;stochastic=true)
    compress(C,D,:L;stochastic=true)

end

function asym_operator_2bond(C,pos,Op,Op_...)
  L=size(C)[1]
  L_uneven = ceil(Int64,(L-1)/2) #amount of uneven bonds
  L_even = floor(Int64,(L-1)/2)  #amout of even bonds
  if pos==:uneven
    (C[1],C[2])=operator_2bond(C[1],C[2],Op_[1]);
    if L%2==0
      (C[2*L_uneven-1],C[2*L_uneven])=operator_2bond(C[2*L_uneven-1],C[2*L_uneven],Op_[2]);
      for j=2:L_uneven-1
        (C[2*j-1],C[2*j])=operator_2bond(C[2*j-1],C[2*j],Op);
      end
    else
      for j=2:L_uneven
        (C[2*j-1],C[2*j])=operator_2bond(C[2*j-1],C[2*j],Op);
      end
    end

  elseif pos==:even
    if L%2==1
      (C[2*L_even],C[2*L_even+1])=operator_2bond(C[2*L_even],C[2*L_even+1],Op_[1]);
      for j=1:L_even-1
        (C[2*j],C[2*j+1])=operator_2bond(C[2*j],C[2*j+1],Op);
      end
    else
      for j=1:L_even
        (C[2*j],C[2*j+1])=operator_2bond(C[2*j],C[2*j+1],Op);
      end
    end
  end
end
