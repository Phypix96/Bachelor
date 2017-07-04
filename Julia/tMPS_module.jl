module MPS

export tMPS, operator, overlap, rand_MPS, compress

"""
  **compress(C, D::Int64=Inf; direction=:R)**

**C**   of Type Array{Any}, every element contains the set of d Matrices for one site\\
**D**   maximum dimension of the Matrices in C, bigger oes will be truncated by a SVD\\
**direction** determines if the Matrices are compressed from left to right (="L") or from
right to left (=:R)\\
  \\
The compressed Array C is normalised, as after the last SVD the scalar U*S (or S*V†
respectively) are ignored.\\
\\
FUNCTION: *Changes the entries of C to normalised, truncated matrices*
Starting at one end of the MPS, the first matrix is subjected to a SVD. Considering
direction=:L, the new matrix is U, while the next matrix in line is multiplied by
SV† and then in turn subjected to a SVD and so on. For the last matrix, SV† is a scalar, that
is ignored to maintain normalisation. If the columnlength of U or V is greater than D,
the matrices are truncated by shrinking the matrices to the first D columns.
"""
function compress(C, D::Int64=Inf, direction=:R)
  L=size(C)[1];
  if direction==:L
    for i=1:L-1
      (C[i],C[i+1])=left_norm(C[i], C[i+1], D);
    end
    C[L]=left_norm(C[L]) #last step normalises C

  elseif direction==:R
    for i=0:L-2
      (C[L-i],C[L-i-1])=right_norm(C[L-i], C[L-i-1], D);
    end
  C[1]=right_norm(C[1]) #last step normalises C
  end
end

"""
  **left_norm(A::Array{Complex128}, M::Array{Complex128}, D::Int64=Inf)**

**A**   matrix that will be left-normalised and truncated\\
**M**   matrix right to A\\
**D**   maximum dimension of A after truncation\\
\\
FUNCTION: *returns two matrices to replace A and M*\\
Set of A-matrices are arranged as a column vector and then subjected to a SVD.
If the size of S is bigger than D, only the first D columns of U and V will be used.
Truncated U and SV†M will be returned.
"""
function left_norm(A::Array{Complex128}, M::Array{Complex128}, D::Int64=Inf)
  (D1,D2,d)=size(A);
  (D3,D4,d)=size(M);
  H=complex(zeros(D1*d, D2));
  #grouping σ-matricies atop each other for svd
  for i = 1:D1 j=1:D2
    for σ =1:d
      H[σ*D1+i-D1,j]=A[i,j,σ];
    end
  end
  (U,S,V)=svd(H);

  D_is=Int64(min(D,size(S)[1]))  #if size has to be truncated, match the matrix dimensions
  M2=complex(zeros(D_is,D4,d))

  for σ = 1:d
    for i=1:D1 j = 1:D_is
      A[i,j,σ]=U[σ*D1+i-D1,j]; #splitting U-maticies at [σ*D1+i-D1,j] into A-matricies again
    end
    for i=1:D3, j=1:D4
      for k=1:D_is
        M2[k,j,σ]+=S[k]*conj(V[i,k])*M[i,j,σ]; #new matrix left of the compressed
      end
    end
  end

  return (A[:,1:D_is,:],M2)
end

"""
  **left_norm(A::Array{Complex128})**

**A**   matrices of right-most site\\
\\
FUNCTION: *returns a normalised matrix to replace A*\\
Subjecting A to a SVD, only U is kept, as SV† is the norm of ΣA_σ' A_σ
"""
function left_norm(A::Array{Complex128})
  (D1,D2,d)=size(A);
  H=complex(zeros(D1*d, D2));
  #grouping σ-matricies atop each other for svd
  for i = 1:D1
    for σ =1:d
      H[σ*D1+i-D1,:]=A[i,:,σ];
    end
  end
  (U,S,V)=svd(H);

  for σ = 1:d
    for i=1:D1
      A[i,1,σ]=U[σ*D1+i-D1,1]; #splitting V^†-maticies at [i,σ*D2+j-D2] into A-matricies again
    end
  end
  return (A)
end


"""
  **right_norm(B::Array{Complex128}, M::Array{Complex128}, D::Int64=Inf)**

**B**   matrix that will be left-normalised and truncated\\
**M**   matrix left to B\\
**D**   maximum dimension of B after truncation\\
\\
FUNCTION: *returns two matrices to replace B and M*\\
Set of B-matrices are arranged as a row vector and then subjected to a SVD.
If the size of S is bigger than D, only the first D columns of U and V will be used.
Truncated MUS and V† will be returned.
"""
function right_norm(B::Array{Complex128}, M::Array{Complex128}, D::Int64=Inf)
  (D1,D2,d)=size(B);
  (D3,D4,d)=size(M)
  H=complex(zeros(D1, D2*d));  #grouping σ-matricies next to each other for svd
  for j = 1:D2 i=1:D1
    for σ =1:d
      H[i,σ*D2+j-D2]=B[i,j,σ];
    end
  end
  (U,S,V)=svd(H);

  D_is=Int64(min(D,size(S)[1]))  #if size has to be truncated, match the matrix dimensions
  M2=complex(zeros(D3,D_is,d))
  for σ = 1:d
    for i=1:D_is j = 1:D2
      B[i,j,σ] = conj(V[σ*D2+j-D2,i]); #splitting V^†-maticies at [i,σ*D2+j-D2] into B-matricies again
    end
    for i=1:D3, j=1:D4
      for k=1:D_is
        M2[i,k,σ] += M[i,j,σ]*U[j,k]*S[k]; #new matrix left of the compressed
      end
    end
  end
  return (B[1:D_is,:,:],M2)
end

"""
  **right_norm(B::Array{Complex128})**

**B**   matrices of left-most site\\
\\
FUNCTION: *returns a normalised matrix to replace B*\\
Subjecting B to a SVD, only V† is kept, as US is the norm of ΣB^σ*B^σ'
"""
function right_norm(B::Array{Complex128})
  (D1,D2,d)=size(B);
  H=complex(zeros(D1, D2*d));  #grouping σ-matricies next to each other for svd
  for j = 1:D2 σ =1:d
    H[:,σ*D2+j-D2]=B[:,j,σ];
  end
  Vt=svd(H)[3]';

  for σ = 1:d
    for j = 1:D2
      B[1,j,σ]=Vt[1,σ*D2+j-D2]; #splitting V^†-maticies at [i,σ*D2+j-D2] into B-matricies again
    end
  end
  return (B)
end

"""
  **tMPS(C, h, t::Float64, Δt::Float64, D::Int64=Inf, dom=:real)**

**C**   of Type Array{Any}, every element contains the set of d matrices for one site\\
**h**   loacal Hamiltonian describing next-neighbour interaction\\
**t**   total time after which C is returned (if expectation-values ect. want to be calculated)\\
**Δt**  time for one Trotter-step\\
**D**   maximum dimension of the matrices in C, bigger oes will be truncated by a SVD\\
**dom**   dertermines if real or imaginary tMPS is calculated
\\
WARNING: The matrices in C have to be complex to avoid an InexactError. If Δt>t,
one Trotter-step will be performed with time Δt. D will be set to inf if not specified.

FUNCTION: *changes the entries of C to compressed matrices, that are time-evolved by h
for a time t*\\
A third-order Trotter decomposition (exp(-im*h_odd*Δt/2)exp(-im*h_even*Δt)exp(-im*h_odd*Δt/2)) is
alternatingly applied to the odd and even bonds of the MPS and compressed after every step.
"""
function tMPS(C, h, t::Float64, Δt::Float64, D::Int64=Inf, dom=:real)
    N=div(t,Δt);
    d=size(C[1])[3];

    if dom==:real
      Op_1=expm(-im*h*Δt/2.)
      Op_2=expm(-im*h*Δt)
    elseif dom==:imag
      Op_1=complex(expm(-h*Δt/2.))
      Op_2=complex(expm(-h*Δt))
    end

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
        compress(C,D,:R)
        compress(C,D,:L)
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
    compress(C,D,:R)
    compress(C,D,:L)

end




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
  **overlap(Φ::Array{Any,1},Ψ::Array{Any,1},domain)**

**Φ**   a MPS\\
**Ψ**   a (different) MPS\\
**domain**    determines whether real(), imag() or abs() of the result is returned\\
\\
FUNCTION:*returns value of* ⟨Φ|Ψ⟩
"""
function overlap(Φ::Array{Any,1},Ψ::Array{Any,1},domain)
  overlap=1.;
  for i=1:size(Ψ)[1]
    sum=0.
    for σ=1:size(Ψ[1])[3]
      sum += *(Φ[i][:,:,σ]',overlap,Ψ[i][:,:,σ])
    end
    overlap = sum
  end

  if domain==:real
    return real(overlap)
  elseif domain ==:imag
    return imag(overlap)
  elseif domain==:abs
    return abs(overlap)
  else
    return overlap
  end
end

"""
  **overlap(Φ::Array{Any,1},O::Array{Any,1},Ψ::Array{Any,1},domain)**

**Φ**   a MPS\\
**O**   a MPO\\
**Ψ**   a (different) MPS\\
**domain**    determines whether real(), imag() or abs() of the result is returned\\
\\
FUNCTION:*returns value of* ⟨Φ|O|Ψ⟩
"""
function overlap(Φ::Array{Any,1},O::Array{Any,1},Ψ::Array{Any,1},domain)
  overlap=1.;
  for i=1:size(Ψ)[1]
    sum=0.
    for σ=1:size(Ψ[1])[3]
      sum += *(Φ[i][:,:,σ]',overlap,operator(Ψ[i],O[i])[:,:,σ])
    end
    overlap = sum
  end

  if domain==:real
    return real(overlap)
  elseif domain ==:imag
    return imag(overlap)
  elseif domain==:abs
    return abs(overlap)
  else
    return overlap
  end
  #domain==:real ? return(real(overlap)):domain==:imag?return(imag(overlap)):domain==:abs?return(abs(overlap)):return(overlap)
end

"""
  **overlap(Φ,Ψ,domain,site_Op::Tuple{Array{},Int64}...)**

**Φ**   a MPS\\
**Ψ**   a (different) MPS\\
**domain**    determines whether real(), imag() or abs() of the result is returned\\
**site_OP**   tuple consisting of a MPO and the site it acts on\\
\\
WARNING: Only one operator can be applied per site, if muliple operators are needed,
calculate the MPO product beforehand.\\

FUNCTION:*returns value of* ⟨Φ|O_i1×O_i2×...|Ψ⟩\\
If a site of a MPO is identical to the current contraction-position, the MPO is
applied, otherwise the usual contraction is
"""
function overlap(Φ,Ψ, domain,site_Op::Tuple{Array{},Int64}...)
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
    end

    if !contraction
      for σ=1:size(Ψ[1])[3]
        sum += *(Φ[i][:,:,σ]',overlap,Ψ[i][:,:,σ])
      end
      overlap = sum
    end

  end
  if domain==:real
    return real(overlap)
  elseif domain ==:imag
    return imag(overlap)
  elseif domain==:abs
    return abs(overlap)
  else
    return overlap
  end

end

"""
  **rand_MPS(C::Array{Any},d::Int64,D::Int64)**

**C**   Array, which will be filled with random matrices in the MPS format\\
**d**   determines how many matrices per site are generated\\
**D**   Maximal Dimension of the matrices\\
\\
FUNCTION:*Filles Array C with random Matrices up to dimension D of the form (1xd)(dxd^2)...(d^2xd)(dx1)*
"""
function rand_MPS(C::Array{Any},d::Int64,D::Int64)
  L=size(C)[1]

  for i=1:div(L,2)
      C[i]=rand(Complex128,min(D,2^(i-1)),min(D,2^i),d)
  end
  if L%2 == 1
    C[div(L,2)+1]=rand(Complex128,min(D,2^(div(L,2))),min(D,2^(div(L,2))),d)
  end
  for i=1:div(L,2)
      C[L-(i-1)]=rand(Complex128,min(D,2^i),min(D,2^(i-1)),d)
  end
  compress(C,D)
end

end
