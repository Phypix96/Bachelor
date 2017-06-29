include("compress.jl")
"""
  **tMPS(C, h, t::Float64, Δt::Float64, D::Int64=Inf)**

**C**   of Type Array{Any}, every element contains the set of d matrices for one site\\
**h**   loacal Hamiltonian describing next-neighbour interaction\\
**t**   total time after which C is returned (if expectation-values ect. want to be calculated)\\
**Δt**  time for one Trotter-step\\
**D**   maximum dimension of the matrices in C, bigger oes will be truncated by a SVD\\
\\
WARNING: The matrices in C have to be complex to avoid an InexactError. If Δt>t,
one Trotter-step will be performed with time Δt. D will be set to inf if not specified.

FUNCTION: *changes the entries of C to compressed matrices, that are time-evolved by h
for a time t*\\
A third-order Trotter decomposition (exp(-im*h_odd*Δt/2)exp(-im*h_even*Δt)exp(-im*h_odd*Δt/2)) is
alternatly applied to the odd and even bonds of the MPS and compressed after every step.
"""
function tMPS(C, h, t::Float64, Δt::Float64, D::Int64=Inf)
    N=floor(t/Δt);
    d=size(C[1])[3];
    Op_1=expm(-im*h*Δt/2.)
    Op_2=expm(-im*h*Δt)
    L_uneven = ceil(Int64,(size(C)[1]-1)/2) #amount of uneven bonds
    L_even = floor(Int64,(size(C)[1]-1)/2)  #amout of even bonds

    for j=1:L_uneven
      (C[2*j-1],C[2*j])=operator_2bond(C[2*j-1],C[2*j],Op_1);
    end

    #as Op_1 commutes with itself, if N>1 two consecutive steps with Op_1 can be combined
    for i=1:N-1

        for j=1:L_even
          (C[2*j],C[2*j+1])=operator_2bond(C[2*j],C[2*j+1],Op_2);
        end
        compress(C,D,direction="R")
        compress(C,D,direction="L")
        for j=1:L_uneven
          (C[2*j-1],C[2*j])=operator_2bond(C[2*j-1],C[2*j],Op_2);
        end
    end

    for j=1:L_even
      (C[2*j],C[2*j+1])=operator_2bond(C[2*j],C[2*j+1],Op_2);
    end
    for j=1:L_uneven
      (C[2*j-1],C[2*j])=operator_2bond(C[2*j-1],C[2*j],Op_1);
    end
    compress(C,D,direction="R")
    compress(C,D,direction="L")

end


"""
  **operator_2bond(A::Array{Complex128},B::Array{Complex128}, Op::Array{Complex128})**

**A** left Matrix of the bond\\
**B** right Matrix of the bond\\
**Op**  Operator applied to a specific bond\\
\\
FUNCTION:*returns two matrices transformed by nearest-neighbour interaction Op to replace A and B*\\
The indices of the operator are regrouped and the operator divided into the product of
two matrices via SVD. The resulting matrices are applied as MPOs to the left and right
sites, A and B, respectively.

"""
function operator_2bond(A::Array{Complex128},B::Array{Complex128}, Op::Array{Complex128})
    (D1,D2,d)=size(A)
    (D1_B,D2_B)=size(B)[1:2]
    P=complex(zeros(d^2,d^2))
    N1=complex(zeros(D1,d^2*D2,d))
    N2=complex(zeros(d^2*D1_B,D2_B,d))
    for i=1:d
      for j=1:d
        for l=1:d
          for k=1:d
            P[i*d+j-d, l*d+k-d] = Op[i*d+l-d, j*d+k-d]
          end
        end
      end
    end

    F=svdfact(P);

    U1=complex(zeros(d^2,d^2))
    U2=complex(zeros(d^2,d^2))
    for i=1:d^2, j=1:d^2
      U1[i,j]=F[:U][i,j]*sqrt(F[:S][j])
      U2[i,j]=sqrt(F[:S][i])*F[:Vt][i,j]
    end

    for σ=1:d
      for i=1:D1
        for j=1:D2
          for k=1:d^2
            for σi=1:d
              N1[i,k*D2+j-D2,σ] += U1[σ*d+σi-d,k]*A[i,j,σi];
            end
          end
        end
      end
    end
    for σ=1:d
      for i=1:D1_B
        for j=1:D2_B
          for k=1:d^2
            for σi=1:d
              N2[k*D1_B+i-D1_B,j,σ] += U2[k,σ*d+σi-d]*B[i,j,σi];
            end
          end
        end
      end
    end
    return (N1,N2)
end

"""
  **operator(A,H)**

**A**   set of matrices on which the operator will be applied\\
**O**   operator applied to A\\
\\
FUNCTION:*returns set of matrices corresponding to O applied on A*\\
applies H to A according to MPO-formalism, returning matrices of dimension dim(O)\*dim(A)
"""
function operator(A,O)
  (A1,A2,d)=size(A)
  (B1,B2)=size(O)[1:2]
  N=complex(zeros(B1*A1,B2*A2,d))
  for a1=1:A1
    for b1=1:B1
      for a2=1:A2
        for b2=1:B2
          for σ1=1:d
            for σ2=1:d
              N[b1*A1+a1-A1,b2*A2+a2-A2,σ1]+=O[b1,b2,σ1,σ2]*A[a1,a2,σ2]
            end
          end
        end
      end
    end
  end
  return N
end

"""
  **overlap(Φ::Array{Any,1},Ψ::Array{Any,1})**

**Φ**   a MPS\\
**Ψ**   a (different) MPS\\
\\
FUNCTION:*returns value of* ⟨Φ|Ψ⟩
"""
function overlap(Φ::Array{Any,1},Ψ::Array{Any,1})
  overlap=1.;
  for i=1:size(Ψ)[1]
    sum=0.
    for σ=1:size(Ψ[1])[3]
      sum += *(Φ[i][:,:,σ]',overlap,Ψ[i][:,:,σ])
    end
    overlap = sum
  end
  return overlap
end

"""
  **overlap(Φ::Array{Any,1},O::Array{Any,1},Ψ::Array{Any,1})**

**Φ**   a MPS\\
**O**   a MPO\\
**Ψ**   a (different) MPS\\
\\
FUNCTION:*returns value of* ⟨Φ|O|Ψ⟩
"""
function overlap(Φ::Array{Any,1},O::Array{Any,1},Ψ::Array{Any,1})
  overlap=1.;
  for i=1:size(Ψ)[1]
    sum=0.
    for σ=1:size(Ψ[1])[3]
      sum += *(Φ[i][:,:,σ]',overlap,operator(Ψ[i],O[i])[:,:,σ])
    end
    overlap = sum
  end
  return overlap
end

"""
  **overlap(Φ,Ψ,site_Op::Tuple{Array{},Int64}...)**

**Φ**   a MPS\\
**Ψ**   a (different) MPS\\
**site_OP**   tuple consisting of a MPO and the site it acts on\\
\\
WARNING: Only one operator can be applied per site, if muliple operators are needed,
calculate the MPO product beforehand. CURRENTLY ONLY ONE OPERATOR WORKS!!!!

FUNCTION:*returns value of* ⟨Φ|O_i1×O_i2×...|Ψ⟩
If a site of a MPO is identical to the current contraction-position, the MPO is
applied, otherwise the usual contraction
"""
function overlap(Φ,Ψ,site_Op::Tuple{Array{},Int64}...)
  overlap = 1.;
  for i = 1:size(Ψ)[1]
    sum = 0.
    contraction = false

    for j = 1:length(site_Op)
      if site_Op[j][2] == i
        contraction = true
        for σ=1:size(Ψ[1])[3]
          sum += *(Φ[i][:,:,σ]',overlap,operator(Ψ[i],site_Op[j][1])[:,:,σ])
        end
        overlap = sum
    end

    if !contraction
        for σ=1:size(Ψ[1])[3]
          sum += *(Φ[i][:,:,σ]',overlap,Ψ[i][:,:,σ])
        end
        overlap = sum
      end
    end

  end

  return overlap
end
