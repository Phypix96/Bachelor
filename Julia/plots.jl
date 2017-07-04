using Plots
using DataFrames
pyplot()

data=readtable("Sz_11.data",separator=',',header=false)
data2=readtable("11.data",separator=',',header=false)
A=zeros(11,201)
A2=zeros(11,201)
#Matrizen nur nötig, da die Differenz der beiden Daten berechne wird, was mit DataFrames nicht geht
for i=1:11
  for j=1:201
    A[i,j]=data[j][i]
    A2[i,j]=data2[i+2][j]
  end
end

for site=1:4
  plt=plot(data2[1],A[site,:],label="tMPS (D=50)")
  plot!(xaxis=("Time",font(14, "Helvetica")), yaxis=("⟨S_z⟩",(-0.5,0.5),font(14, "Helvetica")), title=string("⟨S_z⟩ at site ",site),titlefont=font(16,"Helvetica"))
  plot!(data2[1],data2[site+2],label="Reference Data",wsize=(1020,765))
  savefig(string("Site_",site))
end

heatmap(data2[1],[1:11],A, seriescolor=:inferno,xlabel="Time",ylabel="Sites",title="⟨S_z⟩ with tMPS",wsize=(765,567))
savefig("Heatmap")
heatmap(data2[1],[1:11],A-A2, seriescolor=:viridis,xlabel="Time",ylabel="Sites",title="Difference between tMPS and Reference Data",wsize=(765,567))
savefig("Heatmap_diff")
site=7

plt=plot(data2[1],data[7]-data2[9],label="tMPS")
plot!(xaxis=("Time",font(14, "Helvetica")), yaxis=("⟨S_z⟩",font(14, "Helvetica")), title=string("Δ⟨S_z⟩ at site 7 with full precision"),titlefont=font(16,"Helvetica"))
plot!(wsize=(1020,765))
savefig(string("Site_7_diff"))
