"""
  **compress(C, D::Int64=Inf; direction="R")**

**C**   of Type Array{Any}, every element contains the set of d Matrices for one site\\
**D**   maximum dimension of the Matrices in C, bigger oes will be truncated by a SVD\\
**direction** determines if the Matrices are compressed from left to right (="L") or from
right to left (="R")\\
  \\
The compressed Array C is normalised, as after the last SVD the scalar U*S (or S*V†
respectively) are ignored.\\
\\
FUNCTION: *Changes the entries of C to normalised, truncated matrices*
Starting at one end of the MPS, the first matrix is subjected to a SVD. Considering
direction="L", the new matrix is U, while the next matrix in line is multiplied by
SV† and then in turn subjected to a SVD and so on. For the last matrix, SV† is a scalar, that
is ignored to maintain normalisation. If the columnlength of U or V is greater than D,
the matrices are truncated by shrinking the matrices to the first D columns.
"""
function compress(C, D::Int64=Inf; direction="R")
  L=size(C)[1];
  if direction=="L"
    for i=1:L-1
      (C[i],C[i+1])=left_norm(C[i], C[i+1], D);
    end
    C[L]=left_norm(C[L]) #last step normalises C

  elseif direction=="R"
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
        M2[k,j,σ]+=S[k]*V'[k,i]*M[i,j,σ]; #new matrix left of the compressed
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
  F=svdfact(H);

  for σ = 1:d
    for j = 1:D2
      B[1,j,σ]=F[:Vt][1,σ*D2+j-D2]; #splitting V^†-maticies at [i,σ*D2+j-D2] into B-matricies again
    end
  end
  return (B)
end
