include("time_evolve.jl")
using Plots

L=64
dim=100

t=1.
Δt=0.5

J=1
Jz=1
b=0;
h=[-b/2+Jz/4 0 0 0; 0 -Jz/4 J/2 0; 0 J/2 -Jz/4 0; 0 0 0 Jz/4+b/2]
S_z=[0.5 0;0 -0.5]
S_x=1/2*[0 1;1 0]
S_y=1/2im*complex([0 1;-1 0])
Sx=reshape(S_x,1,1,2,2)
Sy=reshape(S_y,1,1,2,2)
Sz=reshape(S_z,1,1,2,2)
energy=zeros(L-1)

D=Array{Any}(L)
rand_MPS(D,2,dim)

err=1.;
i=1.;
@time tMPS(D,h,20.,2.,dim,:imag)
en_center1=(overlap(D,D,:real,(Sz,div(L,2)),(Sz,div(L,2)+1))+overlap(D,D,:real,(Sx,div(L,2)),(Sx,div(L,2)+1))+overlap(D,D,:real,(Sy,div(L,2)),(Sy,div(L,2)+1)))[1]
while err>1e-5
  @time tMPS(D,h,50.*i,5.*i,dim,:imag)
  en_center2=(overlap(D,D,:real,(Sz,div(L,2)),(Sz,div(L,2)+1))+overlap(D,D,:real,(Sx,div(L,2)),(Sx,div(L,2)+1))+overlap(D,D,:real,(Sy,div(L,2)),(Sy,div(L,2)+1)))[1]
  err=abs(en_center2-en_center1);
  en_center1=en_center2;
  i=i/10.;
  println(err)
end
#@time tMPS(D,h,0.5,0.05,dim,:imag)
#@time tMPS(D,h,0.05,0.005,dim,:imag)
#@time tMPS(D,h,0.005,0.0005,dim,:imag)
for j=1:L-1
  energy[j]=(overlap(D,D,:real,(Sz,j),(Sz,j+1))+overlap(D,D,:real,(Sx,j),(Sx,j+1))+overlap(D,D,:real,(Sy,j),(Sy,j+1)))[1]
end

en=0.;
for i=1:L-1
  en+=energy[i]/(L-1)
end
