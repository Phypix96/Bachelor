include("compress.jl")
include("time_evolve.jl")
using Plots
using DataFrames
pyplot()
data=readtable("11.data",separator=',',header=false)
################################################################################
################################################################################
#main
max_dim=80
L=101
points=50
error=zeros(div(max_dim,5))
#for dim=5:5:max_dim

  dim=50

  site=7
  D=Array{Any}(L)
  D2=Array{Any}(L)
  spin=zeros(points,L)
  energy=zeros(points)

  J=1
  Jz=1
  b=0;
  h=[-b/2+Jz/4 0 0 0; 0 -Jz/4 J/2 0; 0 J/2 -Jz/4 0; 0 0 0 Jz/4+b/2]
  A2=reshape(complex([0.;1.]),1,1,2)
  A1=reshape(complex([1.;0.]),1,1,2)
  S_z=[0.5 0;0 -0.5]
  S_p=[0 1;0 0]
  S_m=[0 0;1 0]
  S_p=reshape(S_p,1,1,2,2)
  S_m=reshape(S_m,1,1,2,2)
  Sz=reshape(S_z,1,1,2,2)
  W1=zeros(1,5,2,2)
  W2=zeros(5,1,2,2)
  W=zeros(5,5,2,2)
  W[1,1,:,:]=eye(2,2)
  W[2,1,:,:]=S_p
  W[3,1,:,:]=S_m
  W[4,1,:,:]=S_z
  W[5,2,:,:]=J/2*S_m
  W[5,3,:,:]=J/2*S_p
  W[5,4,:,:]=Jz*S_z
  W[5,5,:,:]=eye(2,2)
  W1[:,:,:,:]=W[5,:,:,:]
  W2[:,:,:,:]=W[:,1,:,:]
  H=Array{Any}(L)
  for i=2:L-1
    H[i]=W
  end
  H[1]=W1
  H[L]=W2

  for i=1:L
    if i%2==1
      D[i]=A2
    else
      D[i]=A1
    end
  end

  @time println(overlap(D,H,D))
