function compress(C, D::Int64=Inf; direction="R")
  L=size(C)[1];
  norm=0.;
  if direction=="L"
    for i=1:L-1
      (A,M)=left_norm(C[i], C[i+1], D);
      C[i]=A;
      C[i+1]=M;
    end
    for σ=1:size(C[1])[3]
      norm+=*(C[L][:,:,σ]',C[L][:,:,σ])
    end
    #C[L]=C[L]/sqrt(abs(norm[1]))

  elseif direction=="R"
    for i=0:L-2
      (B,M)=right_norm(C[L-i], C[L-i-1], D);
      C[L-i]=B;
      C[L-i-1]=M;
    end

    for σ=1:size(C[1])[3]
      norm+=*(C[1][:,:,σ],C[1][:,:,σ]')
    end
    #C[1]=C[1]/sqrt(abs(norm[1]))
  end
end

function left_norm(A::Array{Complex128}, M::Array{Complex128}, D::Int64=Inf)
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
    M[1:D_is,:,σ]=*(diagm(S[1:D_is]),transpose(V[:,1:D_is]),M[:,:,σ]); #new matrix left of the compressed
  end

  A=A[:,1:min(D,D2),:]
  M=M[1:min(D,size(M)[1]),:,:]
  return (A,M)
end

function right_norm(B::Array{Complex128}, M::Array{Complex128}, D::Int64=Inf)
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
    M[:,1:D_is,σ]=*(M[:,:,σ],U[:,1:D_is],diagm(S[1:D_is])); #new matrix left of the compressed
  end
  B=B[1:min(D,D1),:,:]
  M=M[:,1:min(D,size(M)[2]),:]
  return (B,M)
end
