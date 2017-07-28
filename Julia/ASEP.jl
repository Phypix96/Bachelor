include("time_evolve.jl")

L=10
d=32
t=0.1
Δt=0.05
α=0.2
β=0.2
p_left=0.4
p_right=0.5

h1=[1 0 α 0; 0 1-p_right p_left α; 0 p_right 1-p_left-α 0; 0 0 0 1-α]
h=[1 0 0 0; 0 1-p_right p_left 0; 0 p_right 1-p_left 0; 0 0 0 1]
hl=[1-β 0 0 0; 0 1-p_right-β p_left 0; β p_right 1-p_left 0; 0 β 0 1]

H=[h1 h hl];

D=Array{Any}(L)
for i=1:L
  D[i]=reshape(complex([0;1]),(1,1,2))
end
H=[h1,h,hl]
tMPS(D,H,t,Δt,d,:imag)
