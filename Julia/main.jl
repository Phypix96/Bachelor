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
L=11
points=10
error=zeros(div(max_dim,5))
#for dim=5:5:max_dim

  dim=30

  site=7
  D=Array{Any}(L)
  spin_corr=zeros(points,L-1)
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
  Sp=reshape(S_p,1,1,2,2)
  Sm=reshape(S_m,1,1,2,2)
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


  for i=1:points-1
    for j=1:L-1
      spin_corr[i,j]=abs(overlap(D,D,(Sz,j),(Sz,j+1))[1])
    end
    energy[i]=real(overlap(D,H,D)[1])
    @time tMPS(D,h,0.1,0.1,dim)

    println(i,": ",energy[i])
  end

  for j=1:L-1
    spin_corr[points,j]=real(overlap(D,D,(Sz,j),(Sz,j+1))[1])
  end
  energy[points]=real(overlap(D,H,D)[1])
#print(energy)
"""
  for i=1:L
    error[div(dim,5)]+=((spin[:,i]-data[i+2][1:points])'*(spin[:,i]-data[i+2][1:points]))[1]
  end

  l=@layout([[a;b{0.3h}] c{0.45w}])
  xax_1=("Time", 0:5:20)
  yax_1=("⟨S_z⟩", (-0.5,0.5))
  xax_2=("Matrix Dimension D",(5,max_dim))
  yax_2=("Sum of squared Errors",:log,(0.000001,70))
  xax_3=("Time")
  yax_3=("Site")
  fig=plot(layout=l)
  plot!(fig[1],data[1][1:points],[spin[:,7],data[9][1:points]],xaxis=xax_1, yaxis=yax_1, label=[string("tMPS (D=",dim,")") "Reference Data"],title="⟨S_z⟩ at site 7")
  heatmap!(fig[3],data[1][1:points], [1:L],spin',xaxis=xax_3, yaxis=yax_3,seriescolor=:inferno,title="Heatmap of all Sites")

  plot!(fig[2],[5:5:dim],error[1:div(dim,5)],w=2,xaxis=xax_2, yaxis=yax_2,fill=0, α=0.45,title="Sum of squared Errors",legend=false)
  plot!(wsize=(1100,625))
  savefig(string("Data\\Compression_",div(dim,5)))

end
"""
