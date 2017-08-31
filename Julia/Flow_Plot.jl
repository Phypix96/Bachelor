using Plots;
pyplot();
flow=zeros(201,201)
a=linspace(0,1,201)
b=linspace(0,1,201)
for α=0:0.005:1
  for β=0:0.005:1
    if α<0.5 && α<=β
      flow[Int64(round((α)*200))+1, Int64(round((β)*200))+1]=α*(1-α)
    elseif β<0.5 && β<α
      flow[Int64(round((α)*200))+1, Int64(round((β)*200))+1]=β*(1-β)
    elseif α>=0.5 && β>=0.5
      flow[Int64(round((α)*200))+1, Int64(round((β)*200))+1]=0.25
    end
  end
end

x=[0 0; 0.5 0.5]
y=[0.5 0.5; 0.5 1]
z=[0.5 0.5; 1 0.5]
plt=plot(wsize=(500,500),ylabelrotation=90)
plot!(xaxis=("α",font(12),(0,1),0:0.5:1),yaxis=("β",font(12),(0,1),0:0.5:1))
heatmap!(a,b,flow,seriescolor=:inferno,cbar=:none)
plot!(legend=:none)
plot!(x[:,1],x[:,2],w=3,lcolor=:black)
plot!(y[:,1],y[:,2],w=3,lcolor=:black)
plot!(z[:,1],z[:,2],w=3,lcolor=:black)
plt

savefig("ASEP-Phasen.png")
