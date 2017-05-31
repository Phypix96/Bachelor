using Gallium

function compress(C, D::Int64=Inf; direction="R")
  L=size(C)[1];
  if direction=="L"
    for i=1:L-1
      (A,M)=left_norm(C[i], C[i+1], D);
      C[i]=A;
      C[i+1]=M;
    end
    A=left_norm(C[L],D)
    C[L]=A

  elseif direction=="R"
    for i=0:L-2
      (B,M)=right_norm(C[L-i], C[L-i-1], D);
      C[L-i]=B;
      C[L-i-1]=M;
    end
    B=right_norm(C[1],D)
    C[1]=B
  end
end



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
    #M[k,j,σ]=*(diagm(S[k]),transpose(V[k,i]),M[i,j,σ]);
    for i=1:D3, j=1:D4
      for k=1:D_is
        M2[k,j,σ]+=S[k]*V'[k,i]*M[i,j,σ]; #new matrix left of the compressed
      end
    end
  end

  A=A[:,1:D_is,:]

  return (A,M2)
end

function left_norm(A::Array{Complex128}, D::Int64=Inf)
  (D1,D2,d)=size(A);
  H=complex(zeros(D1*d, D2));
  #grouping σ-matricies atop each other for svd
  for i = 1:D1
    for σ =1:d
      H[σ*D1+i-D1,:]=A[i,:,σ];
    end
  end
  (U,S,V)=svd(H);

  D_is=Int64(min(D,size(S)[1]))  #if size has to be truncated, match the matrix dimensions
  for σ = 1:d
    for i=1:D1 j = 1:D_is
      A[i,j,σ]=U[σ*D1+i-D1,j]; #splitting V^†-maticies at [i,σ*D2+j-D2] into A-matricies again
    end
  end

  A=A[:,1:D_is,:]
  return (A)
end



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
  B2=complex(zeros(D_is,D2,d))
  M2=complex(zeros(D3,D_is,d))
  for σ = 1:d
    for i=1:D_is j = 1:D2
      B2[i,j,σ]=V[σ*D2+j-D2,i]; #splitting V^†-maticies at [i,σ*D2+j-D2] into B-matricies again
    end
    for i=1:D3, j=1:D4
      for k=1:D_is
        M2[i,k,σ]+=M[i,j,σ]*U[j,k]*S[k]; #new matrix left of the compressed
      end
    end
    #[:,1:D_is,σ]=*(M[:,:,σ],U[:,1:D_is],diagm(S[1:D_is])); #new matrix left of the compressed
  end
  return (B2,M2)
end

function right_norm(B::Array{Complex128}, D::Int64=Inf)
  (D1,D2,d)=size(B);
  H=complex(zeros(D1, D2*d));  #grouping σ-matricies next to each other for svd
  for j = 1:D2 σ =1:d
    H[:,σ*D2+j-D2]=B[:,j,σ];
  end
  (U,S,V)=svd(H);

  D_is=Int64(min(D,size(S)[1]))  #if size has to be truncated, match the matrix dimensions
  for σ = 1:d
    for i=1:D_is j = 1:D2
      B[i,j,σ]=V[σ*D2+j-D2,i]; #splitting V^†-maticies at [i,σ*D2+j-D2] into B-matricies again
    end
  end
  B=B[1:D_is,:,:]

  return (B)
end
#breakpoint(left_norm)
