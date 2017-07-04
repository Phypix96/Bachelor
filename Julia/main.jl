include("compress.jl")
include("time_evolve.jl")
using Plots
using DataFrames
pyplot(wsize=(1100,625))
data=readtable("11.data",separator=',',header=false)
################################################################################
################################################################################
#main
min_dim=80
max_dim=80
L=11
points=20
t=0.1
Δt_max=0.1
Δt_min=0.05
site=7
D=Array{Any}(L)
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

error=zeros(div(max_dim,5))
en=plot(xaxis=("Time"),yaxis=("Energy"))

for Δt in (Δt_min,Δt_max)
  for dim=min_dim:5:max_dim
  #dim=60

    for i=1:L
      if i%2==1
        D[i]=A2
      else
        D[i]=A1
      end
    end

    for i=1:points-1
      for j=1:L
        spin[i,j]=overlap(D,D,:real,(Sz,j))[1]
      end
      energy[i]=overlap(D,H,D,:real)[1]
      tMPS(D,h,t,Δt,dim,:real)
      #println(i,": ",energy[i])
    end

    for j=1:L
      spin[points,j]=overlap(D,D,:real,(Sz,j))[1]
    end
    energy[points]=overlap(D,H,D,:real)[1]
    #println(points,": ",energy[points])

    for i=1:L
      error[div(dim,5)]+=((spin[:,i]-data[i+2][1:points])'*(spin[:,i]-data[i+2][1:points]))[1]
    end

    if dim==60 || dim==max_dim
      plot!(en,data[1][1:min(points,101)],energy[1:min(points,101)],label=string("D=",dim," Δt=",Δt))
    end

    l=@layout([[a;b{0.3h}] c{0.45w}])
    xax_1=("Time", 0:5:20)
    yax_1=("⟨S_z⟩", (-0.5,0.5))
    xax_2=("Matrix Dimension D",(5,max_dim))
    yax_2=("Sum of squared Errors",:log,(1e-6,70))
    xax_3=("Time")
    yax_3=("Site")
    fig=plot(layout=l)
    plot!(fig[1],data[1][1:points],[spin[:,7],data[9][1:points]],xaxis=xax_1, yaxis=yax_1, label=[string("tMPS (D=",dim,")") "Reference Data"],title="⟨S_z⟩ at site 7")
    plot!(fig[2],[5:5:dim],error[1:div(dim,5)],w=2,xaxis=xax_2, yaxis=yax_2,fill=0, α=0.45,title="Sum of squared Errors",legend=false)
    heatmap!(fig[3],data[1][1:points], [1:L],spin',xaxis=xax_3, yaxis=yax_3,seriescolor=:inferno,title="Heatmap of all Sites")
    savefig(string("Data\\tMPS_",Δt,"_",dim))
  end
end
