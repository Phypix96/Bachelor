include("compress.jl")
#using Plots
################################################################################
################################################################################
#tMPS

function tMPS(C,h,t,Δt::Float64,D::Int64=Inf)
    N=floor(t/Δt);
    println(N)
    d=size(C[1])[3];
    Op_1=expm(-im*h*Δt/2.)
    Op_2=expm(-im*h*Δt)
    """
    C=tMPS_step(C,Op_1,"odd")
    for i=1:N-1
      C=tMPS_step(C,Op_2,"even")
      C=tMPS_step(C,Op_2,"odd")
      compress(C,D,direction="R")
      compress(C,D,direction="L")
    end
    C=tMPS_step(C,Op_2,"even")
    C=tMPS_step(C,Op_1,"odd")


    for j=1:div(size(C)[1]-1,2)
      (A,B)=operator_2bond(C[2*j-1],C[2*j],Op_1);
      C[2*j-1]=A
      C[2*j]=B
    end
"""
    for i=1:N-1
        for j=1:div(size(C)[1]-1,2)
          (A,B)=operator_2bond(C[2*j],C[2*j+1],Op_2);
          C[2*j]=A
          C[2*j+1]=B
        end
        for j=1:div(size(C)[1]-1,2)
          (A,B)=operator_2bond(C[2*j-1],C[2*j],Op_2);
          C[2*j-1]=A
          C[2*j]=B
        end
        compress(C,D,direction="R")
        compress(C,D,direction="L")
    end

    for j=1:div(size(C)[1]-1,2)
      (A,B)=operator_2bond(C[2*j],C[2*j+1],Op_2);
      C[2*j]=A
      C[2*j+1]=B
    end
    """
    for j=1:div(size(C)[1]-1,2)
      (A,B)=operator_2bond(C[2*j-1],C[2*j],Op_1);
      C[2*j-1]=A
      C[2*j]=B
    end
    """
    compress(C,D,direction="R")
    compress(C,D,direction="L")

end

function operator_en(C,h,D::Int64=Inf)
  C2=Array{Any}(size(C)[1])
  d=size(C[1])[3];
  Op=complex(h)
  for j=1:size(C)[1]-1
    (A,B)=operator_2bond(C[j],C[j+1],Op);
    C2[j]=A
    C2[j+1]=B
  end
  compress(C2,D,direction="R")
  compress(C2,D,direction="L")
  return C2
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

function operator(A,H)
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

  L=11
  dim=50
  J=1
  Jz=1
  b=0;
  h=[-b/2+Jz/4 0 0 0; 0 -Jz/4 J/4 0; 0 J/4 -Jz/4 0; 0 0 0 Jz/4+b/2]
  H=reshape(h,1,1,4,4)
  Sz=[0.5 0;0 -0.5]
  S_p=[0 1;0 0]
  S_m=[0 0;1 0]
  Sz=reshape(Sz,1,1,2,2)
  S_p=reshape(S_p,1,1,2,2)
  S_m=reshape(S_m,1,1,2,2)
  A2=reshape(complex([0.;1.]),1,1,2)
  A1=reshape(complex([1.;0.]),1,1,2)
  spin=zeros(51)
  energ=zeros(51)
  site=1
  D=Array{Any}(L)
  D2=Array{Any}(L)
  W1=zeros(1,5,2,2)
  W2=zeros(5,1,2,2)
  W=zeros(5,5,2,2)
  W[1,1,:,:]=eye(2,2)
  W[2,1,:,:]=S_p
  W[3,1,:,:]=S_m
  W[4,1,:,:]=Sz
  W[5,2,:,:]=J/2*S_m
  W[5,3,:,:]=J/2*S_p
  W[5,4,:,:]=Jz*Sz
  W[5,5,:,:]=eye(2,2)
  W1[:,:,:,:]=W[:,5,:,:]
  W2[:,:,:,:]=W[1,:,:,:]
  H=Any[W1,W,W,W,W,W,W,W,W,W,W2]

  for i=1:L
    if i%2==1
      D[i]=A2
    else
      D[i]=A1
    end
  end

for i=1:5
  for j=1:L
    D2[j]=D[j]
  end
  D2[site]=operator(D[site],Sz)
  #compress(D2,dim,direction="R")
  #compress(D2,dim,direction="L")
  spin[i]=real((overlap(D2,D)/overlap(D,D))[1])
  @time tMPS(D,h,0.2,0.1,dim)

  #for j=1:L
  #  D2[j]=operator(D[j],H[j])
  #end
  #D2=operator_en(D,h,20)
  #V=operator(D[6],Sz)
  #compress(D2,50,direction="L")
  #compress(D2,50,direction="R")
  #D2[6]=V
  #energ[i]=real((overlap(D2,D)/overlap(D,D))[1])

  println(i,": ",spin[i])
end

for j=1:L
  D2[j]=D[j]
end
D2[site]=operator(D[site],Sz)
#compress(D2,dim,direction="R")
#compress(D2,dim,direction="L")
spin[51]=real((overlap(D2,D)/overlap(D,D))[1])
print(spin)
x=linspace(0,3,51)
#plot(x,spin)

#Site 1
#t=1,Δt=0.1
spin1=[-0.5,-0.499346,-0.497383,-0.494117,-0.489554,-0.483706,-0.476586,-0.468214,-0.45861,-0.447799,-0.435811,-0.422677,-0.408431,-0.393113,-0.376764,-0.359427,-0.34115,-0.321982,-0.301976,-0.281185,-0.259667,-0.23748,-0.214684,-0.191342,-0.167516,-0.143272,-0.118674,-0.0937888,-0.0686838,-0.0434258,-0.0180826,0.00727851,0.03259,0.057785,0.0827969,0.10756,0.132009,0.15608,0.179711,0.202839,0.225404,0.247348,0.268615,0.289147,0.308893,0.327801,0.345822,0.362908,0.379016,0.394102,0.408128]
#t=0.4,Δt=0.2
spin2=[-0.5,-0.497508,-0.490057,-0.47772,-0.46062,-0.438928,-0.41286,-0.382677,-0.348682,-0.311213,-0.270646,-0.227388,-0.181871,-0.13455,-0.0858997,-0.0364055,0.0134375,0.0631315,0.11218,0.160094,0.206394,0.250619,0.292329,0.331106,0.366564,0.398348,0.426143,0.44967,0.468694,0.483024,0.492517,0.497077,0.496658,0.491263,0.480946,0.465808,0.446002,0.421725,0.393219,0.360769,0.324701,0.285374,0.243184,0.198552,0.151924,0.103768,0.0545645,0.00480601,-0.0450104,-0.0943873,-0.142832]
#t=1,Δt=0.5
spin3=[-0.5,-0.484639,-0.439435,-0.367171,-0.272433,-0.161272,-0.0407393,0.0816258,0.198292,0.302163,0.386933,0.447396,0.479734,0.481783,0.453252,0.395822,0.3131,0.210369,0.0941915,-0.0280958,-0.148894,-0.260813,-0.357066,-0.431802,-0.480396,-0.499728,-0.488434,-0.447072,-0.378173,-0.286098,-0.176727,-0.0570102,0.0655499,0.183395,0.289347,0.376979,0.440922,0.477156,0.483285,0.458759,0.404998,0.325362,0.224933,0.110129,-0.011784,-0.133214,-0.246719,-0.345412,-0.4233,-0.475579,-0.498914]
#t=0.03,Δt=0.015
spin4=[-0.5,-0.499986,-0.499944,-0.499873,-0.499775,-0.499648,-0.499494,-0.499311,-0.4991,-0.498861,-0.498594,-0.498299,-0.497976,-0.497625,-0.497246,-0.496839,-0.496404,-0.495941,-0.495451,-0.494932,-0.494386,-0.493811,-0.493209,-0.492579,-0.491922,-0.491237,-0.490524,-0.489784,-0.489016,-0.48822,-0.487397,-0.486547,-0.485669,-0.484764,-0.483832,-0.482872,-0.481886,-0.480872,-0.479831,-0.478763,-0.477669,-0.476547,-0.475399,-0.474223,-0.473022,-0.471793,-0.470538,-0.469257,-0.467949,-0.466615,-0.465254,-0.463868,-0.462455,-0.461016,-0.459552,-0.458061,-0.456545,-0.455003,-0.453436,-0.451843,-0.450224,-0.448581,-0.446912,-0.445217,-0.443498,-0.441754,-0.439985,-0.438191,-0.436373,-0.43453,-0.432663,-0.430771,-0.428855,-0.426915,-0.424951,-0.422963,-0.420952,-0.418916,-0.416857,-0.414775,-0.412669,-0.41054,-0.408388,-0.406213,-0.404015,-0.401794,-0.399551,-0.397286,-0.394998,-0.392687,-0.390355,-0.388001,-0.385625,-0.383227,-0.380807,-0.378367,-0.375905,-0.373421,-0.370917,-0.368392,-0.365846,-0.36328,-0.360693,-0.358086,-0.355459,-0.352812,-0.350144,-0.347458,-0.344751,-0.342026,-0.339281,-0.336517,-0.333734,-0.330932,-0.328112,-0.325273,-0.322416,-0.319541,-0.316647,-0.313736,-0.310808,-0.307861,-0.304898,-0.301917,-0.29892,-0.295905,-0.292874,-0.289827,-0.286763,-0.283683,-0.280587,-0.277475,-0.274348,-0.271205,-0.268047,-0.264874,-0.261686,-0.258484,-0.255266,-0.252035,-0.248789,-0.245529,-0.242256,-0.238968,-0.235668,-0.232354,-0.229027,-0.225687,-0.222335,-0.21897,-0.215592,-0.212203,-0.208801,-0.205388,-0.201963,-0.198527,-0.19508,-0.191622,-0.188153,-0.184673,-0.181183,-0.177683,-0.174173,-0.170653,-0.167124,-0.163585,-0.160036,-0.156479,-0.152913,-0.149339,-0.145756,-0.142165,-0.138565,-0.134959,-0.131344,-0.127722,-0.124093,-0.120457,-0.116814,-0.113165,-0.109509,-0.105847,-0.102179,-0.0985053,-0.0948261,-0.0911415,-0.0874519,-0.0837573,-0.080058,-0.0763542,-0.0726461,-0.068934,-0.0652179,-0.0614982,-0.057775,-0.0540486,-0.0503191,-0.0465868,-0.0428519,-0.0391146,-0.0353751,-0.0316336,-0.0278903,-0.0241454,-0.0203992,-0.0166518,-0.0129035,-0.0091545,-0.00540496,-0.00165512,0.00209482,0.00584464,0.00959413,0.0133431,0.0170913,0.0208385,0.0245846,0.0283293,0.0320724,0.0358136,0.0395529,0.04329,0.0470246,0.0507565,0.0544856,0.0582117,0.0619344,0.0656537,0.0693693,0.073081,0.0767886,0.0804918,0.0841906,0.0878846,0.0915736,0.0952575,0.098936,0.102609,0.106276,0.109937,0.113592,0.117241,0.120883,0.124518,0.128147,0.131768,0.135381,0.138987,0.142586,0.146176,0.149758,0.153331,0.156896,0.160452,0.163999,0.167537,0.171065,0.174584,0.178093,0.181592,0.185081,0.188559,0.192027,0.195484,0.19893,0.202364,0.205788,0.209199,0.212599,0.215987,0.219363,0.222727,0.226078,0.229416,0.232742,0.236054,0.239353,0.242638,0.24591,0.249168,0.252413,0.255642,0.258858,0.262059,0.265245,0.268416,0.271572,0.274713,0.277839,0.280948,0.284042,0.28712,0.290182,0.293228,0.296257,0.299269,0.302265,0.305243,0.308205,0.311149,0.314075,0.316984,0.319875,0.322748,0.325603,0.32844,0.331258,0.334057,0.336838,0.3396,0.342342,0.345066,0.34777,0.350454,0.353119,0.355764,0.358388,0.360993,0.363577,0.366141,0.368685,0.371207,0.373709,0.37619,0.378649,0.381087,0.383504,0.385899,0.388273,0.390624,0.392954,0.395262,0.397547,0.39981,0.40205,0.404268,0.406463,0.408636,0.410785,0.412911,0.415014,0.417094,0.41915,0.421182,0.423191,0.425176,0.427137,0.429074,0.430987,0.432876,0.43474,0.43658,0.438396,0.440186,0.441952,0.443694,0.44541,0.447101,0.448767,0.450407,0.452023,0.453613,0.455177,0.456716,0.458229,0.459716,0.461178,0.462613,0.464023,0.465406,0.466763,0.468094,0.469399,0.470677,0.471929,0.473154,0.474353,0.475525,0.47667,0.477788,0.47888,0.479944,0.480982,0.481992,0.482976,0.483932,0.484861,0.485763,0.486637,0.487484,0.488304,0.489096,0.489861,0.490598,0.491307,0.491989,0.492643,0.49327,0.493868,0.494439,0.494982,0.495498,0.495985,0.496444,0.496876,0.49728,0.497655,0.498003,0.498323,0.498614,0.498878,0.499113,0.499321,0.4995,0.499651,0.499774,0.499869,0.499936,0.499975,0.499986,0.499968,0.499923,0.499849,0.499747,0.499617,0.499459,0.499273,0.499059,0.498817,0.498546,0.498248,0.497921,0.497567,0.497185,0.496774,0.496336,0.49587,0.495375,0.494853,0.494304,0.493726,0.49312,0.492487,0.491826,0.491138,0.490422,0.489678,0.488907,0.488108,0.487282,0.486428,0.485547,0.484639,0.483703,0.48274,0.48175,0.480733,0.479689,0.478618,0.47752,0.476395,0.475244,0.474065,0.47286,0.471629,0.47037,0.469086,0.467775,0.466437,0.465074,0.463684,0.462268,0.460827,0.459359,0.457865,0.456346,0.454801,0.45323,0.451634,0.450013,0.448366,0.446694,0.444997,0.443275,0.441527,0.439755,0.437959,0.436137,0.434292,0.432421,0.430527,0.428608,0.426665,0.424698,0.422707,0.420693,0.418654,0.416592,0.414507,0.412399,0.410267]

#Site 6
#t=1,Δt=0.1
spin1=[0.5,0.499099,0.496403,0.491931,0.485716,0.477801,0.468241,0.457098,0.444447,0.430368,0.414948,0.398278,0.380454,0.361574,0.341737,0.321042,0.299586,0.277465,0.25477,0.231589,0.208004,0.184095,0.159935,0.13559,0.111122,0.0865903,0.0620456,0.0375368,0.0131083,-0.0111981,-0.0353433,-0.0592893,-0.0829994,-0.106437,-0.129564,-0.152343,-0.174733,-0.196693,-0.218178,-0.23914,-0.259531,-0.279297,-0.298385,-0.316736,-0.334293,-0.350996,-0.366784,-0.381598,-0.395377,-0.408064,-0.419603]
#t=0.4,Δt=0.2
spin2=[0.5,0.497416,0.489695,0.476935,0.45929,0.436974,0.410249,0.379427,0.344855,0.306917,0.266021,0.222597,0.17709,0.129954,0.0816519,0.0326474,-0.0165943,-0.0656103,-0.113941,-0.161134,-0.206744,-0.250336,-0.29149,-0.3298,-0.364884,-0.396385,-0.423976,-0.447365,-0.466302,-0.480582,-0.490049,-0.494598,-0.494181,-0.488802,-0.478525,-0.463463,-0.443782,-0.419698,-0.391467,-0.359387,-0.323786,-0.285026,-0.243489,-0.199579,-0.153717,-0.106337,-0.0578819,-0.00880523,0.040434,0.0893721,0.137543]
#t=1,Δt=0.5
spin3=[0.5,0.481293,0.42788,0.346634,0.246382,0.135686,0.0214696,-0.0911013,-0.197739,-0.293913,-0.374262,-0.432955,-0.464873,-0.467004,-0.439357,-0.384915,-0.308633,-0.215997,-0.111889,-0.00036768,0.114563,0.227689,0.331721,0.417565,0.475926,0.499665,0.485884,0.436751,0.35871,0.260463,0.150713,0.0366344,-0.0764087,-0.18408,-0.281933,-0.364717,-0.42661,-0.462314,-0.468468,-0.44466,-0.393518,-0.319817,-0.229051,-0.126218,-0.0154612,0.0992782,0.213015,0.318749,0.407561,0.470052,0.498662]
