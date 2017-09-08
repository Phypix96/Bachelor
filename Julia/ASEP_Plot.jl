#mean_density=((density[:,:,1]+density[:,:,2]+density[:,:,3]+density[:,:,4]+density[:,:,5]+density[:,:,6]+density[:,:,7]+density[:,:,8]+density[:,:,9]+density[:,:,10])/10)

font = Plots.font("Serif", 8)
pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font,titlefont=font)
Err=zeros(41,41)
mean_flow=zeros(41,41)
std_flow=zeros(41,41)
mean_density=zeros(41,41)
min_entropy=zeros(41,41)

for i=1:41
  for j=1:41
    for l=1:L
      Err[i,j]+=abs(density[i,j,l]-density_prof[i,j,l])^2/L
    end
    mean_flow[i,j]=mean(flow[i,j,:])
    std_flow[i,j]=std(flow[i,j,1:L-1])/mean_flow[i,j]
    mean_density=mean(density[i,j,:])
    Err[i,j]=sqrt(Err[i,j])
    min_entropy[i,j]=minimum(ent_entropy_MPS[i,j,1:L-1])
  end
end
x1=linspace(0,1,40)
x2=linspace(0.1,1,37)
#contour(x1,x1,flow_prof[2:41,2:41],fill=true,nlevels=50,seriescolor=:magma)
#contour(x1,x1,-flow_prof[2:41,2:41]+mean_flow[2:41,2:41],fill=true,nlevels=50,seriescolor=:magma_r)
#contour(x2,x2,std_flow[5:41,5:41],fill=true,nlevels=50,seriescolor=:magma_r)
#contour(x1,x1,min_entropy[2:41,2:41],fill=true,nlevels=150,seriescolor=cgrad(:magma_r, scale=:exp))
#contour(x1,x1,Err[2:41,2:41],fill=true,nlevels=50,seriescolor=:magma_r)

#contour(x1,x1,imfilter(min_entropy[2:41,2:41],Kernel.gaussian(0.5)),seriescolor=:magma_r,xlims=(0,1),ylims=(0,1),fill=true, nlevels=50)
#contour(x1,x1,imfilter(-flow_prof[2:41,2:41]+mean_flow[2:41,2:41],Kernel.gaussian(0.5)),seriescolor=:magma_r,xlims=(0,1),ylims=(0,1),fill=true, nlevels=50)
contour(x2,x2,imfilter(std_flow[5:41,5:41],Kernel.gaussian(0.5)),seriescolor=:magma_r,xlims=(0.1,1),ylims=(0.1,1),fill=true, nlevels=50)
#contour(x1,x1,imfilter(Err[2:41,2:41],Kernel.gaussian(0.5)),seriescolor=:magma_r,xlims=(0,1),ylims=(0,1),fill=true, nlevels=50)

plot!(xaxis=(" α",0:0.25:1),yaxis=(" β",0:0.25:1),colorbar_title=" ",wsize=(700,500),dpi=150)
#savefig("Data/entropy_exp_scale.png")

#a=zeros(1000); b=zeros(1000);
#for i=1:1000
#  A=rand(2,2)
#  A2=nnmf(A,2)
#  a[i]=maximum(abs(A-*(A2.W,A2.H)))
#end
