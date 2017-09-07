include("time_evolve.jl")
using Plots, LaTeXStrings

L=6
d=8
t=10.
Δt=0.3
iter=3

α=0.2
β=0.7
p_left=0.0
p_right=1.0

N=reshape([1. 0.; 0. 0.],(1,1,2,2))
Ident=reshape([1. 0.; 0. 1.],(1,1,2,2))
N=complex(N); Ident=complex(Ident)
Right=Array{Any}(L)
Left=Array{Any}(L)

ent_entropy_MPS=zeros(41,41,L-1)

density=zeros(41,41,L)
flow=zeros(41,41,L-1)
density_prof=zeros(41,41,L)
flow_prof=zeros(41,41)


for k=1:41
  for l=1:41
α=(k-1)/40
β=(l-1)/40
h1=[0 0 α 0; 0 -p_right p_left α; 0 p_right -p_left-α 0; 0 0 0 -α]
h=[0 0 0 0; 0 -p_right p_left 0; 0 p_right -p_left 0; 0 0 0 0]
hl=[-β 0 0 0; β -p_right p_left 0; 0 p_right -p_left-β 0; 0 0 β 0]

#(h1,h,hl)=(complex(h1),complex(h),complex(hl))
H=[h1,h,hl];
Ht=[h1',h',hl']


#density=zeros(iter+1,L)
#flow=zeros(iter+1,L-1)
#density_prof=zeros(L)

for i=1:L
  #Right[i]=complex(reshape([1/sqrt(2);1/sqrt(2)],(1,1,2))); Left[i]=complex(reshape([1/sqrt(2);1/sqrt(2)],(1,1,2)))
  Right[i]=reshape([0.5;0.5],(1,1,2))
  #Left[i]=reshape([0.0;1.0],(1,1,2))
end

for i=1:iter
  print(i,": ")
  tMPS(Right,H,t,Δt,d)
  #tMPS(Left,Ht,t,Δt,d)
  for j=1:L
    #density[i,j]=overlap(Left,Right,:real,(N,j))[1]/overlap(Left,Right,:real)[1]
    #density[i,j]=stoc_expect_val(Right,(N,j))[1]/stoc_expect_val(Right,(Ident,1))[1]
    if j<L
      #flow[i,j]=overlap(Left,Right,:real,(N,j),(Ident-N,j+1))[1]/overlap(Left,Right,:real,(Ident,1))[1]
    end
  end
end

#tMPS(Right,H,float(L)*3,0.1,d)
#tMPS(Left,Ht,float(L)*3,0.1,d)
#tMPS(Right,H,float(L)*0.25,0.025,d)
#tMPS(Left,Ht,float(L)*0.25,0.025,d)
tMPS(Right,H,float(L),0.2,d)
for j=1:L
  #density[k, l]=overlap(Left,Right,:real,(N,j))[1]/overlap(Left,Right,:real,(Ident,1))[1]
  density[k,l,j]=stoc_expect_val(Right,(N,j))[1]/stoc_expect_val(Right,(Ident,1))[1]
  #density[iter+1,j]=stoc_expect_val(Right,(N,j))[1]/stoc_expect_val(Right,(Ident,1))[1]
  density_prof[k, l,j]=ASEP_density(j,L,α,β)
  if j<L
    #flow[k,l,j]=overlap(Left,Right,:real,(N,j),(Ident-N,j+1))[1]/overlap(Left,Right,:real,(Ident,1))[1]
    #ent_entropy_MPS[k,l,:]=entropy(Right)
    flow[k,l,j]=stoc_expect_val(Right,(N,j),(Ident-N,j+1))[1]/stoc_expect_val(Right,(Ident,1))[1]
    ent_entropy_MPS[k,l,:]=stoc_entropy(Right)
  end
end
#(ent_entropy_MPS[Int64(round((α-0.0)*40))+1, Int64(round((β-0.0)*40))+1],max_ent_entropy_MPS[Int64(round((α-0.0)*40))+1, Int64(round((β-0.0)*40))+1])=entropy(Right)

flow_prof[k,l]=R(L-1,α,β)/R(L,α,β)



println(α," ",β)
end
end
#font = Plots.font("Serif", 16)
#pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font,titlefont=font)

#min_den_err=signif(min(minimum(density_prof-density[iter+1,:]),0),1)
#max_den_err=signif(max(maximum(density_prof-density[iter+1,:]),0),1)
#yax2=(L"$Δτ_i$", min_den_err:(max_den_err-min_den_err)/2:max_den_err ,(1.1*min_den_err,1.1*max_den_err))

#plt1=plot(linspace(1,L,L),density[iter+1,:],xlabel="i",yaxis=(L"$τ_i$",(0,1),0:0.25:1), legend=false,fill=(0,0.25))
#plt2=plot(linspace(1,L,L),density_prof-density[iter+1,:],xlabel="i",yaxis=yax2,legend=false,yformatter=:scientific)
#plot(plt1,plt2,wsize=(500,500),layout=@layout([a{0.8h};b]))
#savefig(string("L=",L,",D_max=",d,",α=",α,",β=",β))
