include("compress.jl")
using Plots
################################################################################
################################################################################
#tMPS

function tMPS(C,h,t,Δt::Float64,D::Int64=Inf)
    N=floor(t/Δt);
    d=size(C[1])[3];
    Op_odd=expm(-im*h*Δt/2)
    Op_even=expm(-im*h*Δt)
    for i=1:N
        for j=1:div(size(C)[1]-1,2)
          (A,B)=operator_2bond(C[2*j-1],C[2*j],Op_odd);
          C[2*j-1]=A
          C[2*j]=B
        end
        for j=1:div(size(C)[1]-1,2)
          (A,B)=operator_2bond(C[2*j],C[2*j+1],Op_even);
          C[2*j]=A
          C[2*j+1]=B
        end
        for j=1:div(size(C)[1]-1,2)
          (A,B)=operator_2bond(C[2*j-1],C[2*j],Op_odd);
          C[2*j-1]=A
          C[2*j]=B
        end
        compress(C,D,direction="R")
        compress(C,D,direction="L")
    end
end

function operator_2bond(A::Array{Complex128},B::Array{Complex128}, Op::Array{Complex128})
    (D1,D2,d)=size(A)
    (D1_B,D2_B,d)=size(B)
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
    U1=*(F[:U],diagm(sqrt(F[:S])));
    U2=*(diagm(sqrt(F[:S])),F[:Vt]);
    #println(U1)
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
    #U1=reshape(U1,1,d^2,d,d)
    #U2=reshape(U2,d^2,1,d,d)
    #N12=operator(A,U1)
    #N22=operator(B,U2)
    A=N1;
    B=N2;
    return (A,B)
end

function operator(A::Array{Complex128},H)
  (A1,A2,d)=size(A)
  (B1,B2)=size(H)[1:2]

  N=complex(zeros(B1*A1,B2*A2,d))
  for a1=1:A1
    for b1=1:B1
      for a2=1:A2
        for b2=1:B2
          for σ1=1:d
            for σ2=1:d
              N[b1*A1+a1-A1,b2*A2+a2-A2,σ1]+=H[b1,b2,σ1,σ2]*A[a1,a2,σ2]
            end
          end
        end
      end
    end
  end
  return N
end

function overlap(Φ,Ψ)
  overlap=1.;
  N=size(Ψ)[1]
  for i=1:N
    sum=0.
    for σ=1:size(Ψ[1])[3]
      sum += *(Φ[i][:,:,σ]',overlap,Ψ[i][:,:,σ])
    end
    overlap = sum
  end
  return overlap
end

################################################################################
################################################################################
#main
#for i=1:10
"""
  #A=complex(randn(1,2,2));
  #B=complex(randn(2,4,2));
  #C=complex(randn(4,2,2));
  #D=complex(randn(2,1,2));
  #G=Any[A,B,C,D]
"""
  L=11
  J=1
  Jz=1
  b=0;
  h=[-b/2+Jz/4 0 0 0; 0 -Jz/4 J/4 0; 0 J/4 -Jz/4 0; 0 0 0 Jz/4+b/2]
  Sz=[0.5 0;0 -0.5]
  S_p=[0 1;0 0]
  S_m=[0 0;1 0]
  Sz=reshape(Sz,1,1,2,2)
  S_p=reshape(S_p,1,1,2,2)
  S_m=reshape(S_m,1,1,2,2)
  A2=reshape(complex([0.;1.]),1,1,2)
  A1=reshape(complex([1.;0.]),1,1,2)
  spin=zeros(21)
  D=Array{Any}(L)
  D2=Array{Any}(L)
  for i=1:L
    if i%2==1
      D[i]=A1
    else
      D[i]=A2
    end
  end
for i=1:21
  @time tMPS(D,h,0.1,0.005,30)
  V=operator(D[6],Sz)
  for j=1:L
    D2[j]=D[j]
  end
  D2[6]=V
  spin[i]=real((overlap(D2,D)/overlap(D,D))[1])
  println(i)
end
print(spin)
x=linspace(0,2,21)
plot(x,spin)
gui()
