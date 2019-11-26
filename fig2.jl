using Distributions
using AdaptiveABC
using Distributed
using DataFrames
using CSV
# 1) make funcs consisent 2)

# model definition (priors over parameters not already fixed by prefit)
#poisson model for n different genes with parameters: constant transcription rate, constant decay rate and b of NLM
model_pois_preset(n)=vcat(fill(Uniform(0,1000),n),fill(Uniform(2,10),n),Uniform(0,40))
#error function wrapper for specific model. n as above, gps the prefit NLM parameters and gr_dep a flag for growth rate dependence of transcription (0 or 1)
function f_pois_preset(n,gps,gr_dep,o=3,out="fit")
  func(x,y)=rho_general(x,vcat(repeat(vcat(1,0),outer=n::Int64),repeat(vcat(0,0),outer=n),vcat(transpose(y[1:n]),transpose(fill(0,n)))[:],vcat(transpose(y[n+1:2*n]),transpose(fill(0,n)))[:],gps[1],y[2*n+1],gps[3:end]),o,x->sqrt(sum(x)),out,gr_dep)
end

#model and wrapper for Poisson with size dependent transcription rate
model_pois_v_preset(n)= vcat(fill(Uniform(0,200),n),fill(Uniform(2,10),n),Uniform(0,40))
function f_pois_v_preset(n,gps,gr_dep,o=3,out="fit")
  func(x,y)=rho_general(x,vcat(repeat(vcat(1,0),outer=n::Int64),repeat(vcat(0,0),outer=n),vcat(transpose(fill(0,n)),transpose(y[1:n]))[:],vcat(transpose(y[n+1:2*n]),transpose(fill(0,n)))[:],gps[1],y[2*n+1],gps[3:end]),o,x->sqrt(sum(x)),out,gr_dep)
end

#model and wrapper for bursty with constant rates
model_bursty_preset(n)= vcat(fill(Uniform(0,100),n),fill(Uniform(0,10),n),fill(Uniform(2,10),n),Uniform(0,40))
function f_bursty_preset(n,gps,gr_dep,o=3,out="fit")
  func(x,y)=rho_bursty(x,vcat(vcat(transpose(y[1:n]),transpose(fill(0,n)))[:],vcat(transpose(y[n+1:2*n]),transpose(fill(0,n)))[:],vcat(transpose(y[2*n+1:3*n]),transpose(fill(0,n)))[:],gps[1],y[3*n+1],gps[3:end]),o,x->sqrt(sum(x)),out,gr_dep)
end

#model and wrapper for bursty with size dependent bursty size
model_bursty_v_preset(n)= vcat(fill(Uniform(0,100),n),fill(Uniform(0,2),n),fill(Uniform(2,10),n),Uniform(0,40))
function f_bursty_v_preset(n,gps,gr_dep,o=3,out="fit")
  func(x,y)=rho_bursty(x,vcat(vcat(transpose(y[1:n]),transpose(fill(0,n)))[:],vcat(transpose(fill(0,n)),transpose(y[n+1:2*n]))[:],vcat(transpose(y[2*n+1:3*n]),transpose(fill(0,n)))[:],gps[1],y[3*n+1],gps[3:end]),o,x->sqrt(sum(x)),out,gr_dep)
end

#model and wrapper for bursty with size dependent burst frequency
model_bursty_kon_preset(n)= vcat(fill(Uniform(0,20),n),fill(Uniform(0,10),n),fill(Uniform(2,10),n),Uniform(0,40))
function f_bursty_kon_preset(n,gps,gr_dep,o=3,out="fit")
  func(x,y)=rho_bursty(x,vcat(vcat(transpose(fill(0,n)),transpose(y[1:n]))[:],vcat(transpose(y[n+1:2*n]),transpose(fill(0,n)))[:],vcat(transpose(y[2*n+1:3*n]),transpose(fill(0,n)))[:],gps[1],y[3*n+1],gps[3:end]),o,x->sqrt(sum(x)),out,gr_dep)
end

#model and wrapper for Poisson with cell cycle effects and constant rates
model_pois_cc_preset(n)= vcat(fill(Uniform(0,1000),n),fill(Uniform(2,10),n),fill(Uniform(0,1),3*n),Uniform(0,40))
function f_pois_cc_preset(n,gps,gr_dep,o=3,out="fit")
  func(x,y)=rho_cc(x,vcat(repeat(vcat(1,0),outer=n::Int64),repeat(vcat(0,0),outer=n),vcat(transpose(y[1:n]),transpose(fill(0,n)))[:],vcat(transpose(y[n+1:2*n]),transpose(fill(0,n)))[:],y[2*n+1:2*n+3],gps[1],y[end],gps[3:end]),o,x->sqrt(sum(x)),out,gr_dep)
end

#model and wrapper for Poisson with cell cycle effects and size dependent transctiption rate
model_pois_cc_v_preset(n)= vcat(fill(Uniform(0,200),n),fill(Uniform(2,10),n),fill(Uniform(0,1),3*n),Uniform(0,40))
function f_pois_cc_v_preset(n,gps,gr_dep,o=3,out="fit")
  func(x,y)=rho_cc(x,vcat(repeat(vcat(1,0),outer=n::Int64),repeat(vcat(0,0),outer=n),vcat(transpose(fill(0,n)),transpose(y[1:n]))[:],vcat(transpose(y[n+1:2*n]),transpose(fill(0,n)))[:],y[2*n+1:2*n+3],gps[1],y[end],gps[3:end]),o,x->sqrt(sum(x)),out,gr_dep)
end

#model and wrapper for bursty with cell cycle effects and constant rates
model_bursty_cc_preset(n)= vcat(fill(Uniform(0,100),n),fill(Uniform(0,10),n),fill(Uniform(2,10),n),fill(Uniform(0,1),3*n),Uniform(0,40))
function f_bursty_cc_preset(n,gps,gr_dep,o=3,out="fit")
  func(x,y)=rho_bursty_cc(x,vcat(vcat(transpose(y[1:n]),transpose(fill(0,n)))[:],vcat(transpose(y[n+1:2*n]),transpose(fill(0,n)))[:],vcat(transpose(y[2*n+1:3*n]),transpose(fill(0,n)))[:],y[3*n+1:3*n+3],gps[1],y[end],gps[3:end]),o,x->sqrt(sum(x)),out,gr_dep)
end

#model and wrapper for bursty with cell cycle effects and size dependent burst size
model_bursty_cc_v_preset(n)= vcat(fill(Uniform(0,100),n),fill(Uniform(0,2),n),fill(Uniform(2,10),n),fill(Uniform(0,1),3*n),Uniform(0,40))
function f_bursty_cc_v_preset(n,gps,gr_dep,o=3,out="fit")
  func(x,y)=rho_bursty_cc(x,vcat(vcat(transpose(y[1:n]),transpose(fill(0,n)))[:],vcat(transpose(fill(0,n)),transpose(y[n+1:2*n]))[:],vcat(transpose(y[2*n+1:3*n]),transpose(fill(0,n)))[:],y[3*n+1:3*n+3],gps[1],y[end],gps[3:end]),o,x->sqrt(sum(x)),out,gr_dep)
end

#model and wrapper for bursty with cell cycle effects and size dependent burst frequency
model_bursty_cc_kon_preset(n)= vcat(fill(Uniform(0,20),n),fill(Uniform(0,10),n),fill(Uniform(2,10),n),fill(Uniform(0,1),3*n),Uniform(0,40))
function f_bursty_cc_kon_preset(n,gps,gr_dep,o=3,out="fit")
  func(x,y)=rho_bursty_cc(x,vcat(vcat(transpose(fill(0,n)),transpose(y[1:n]))[:],vcat(transpose(y[n+1:2*n]),transpose(fill(0,n)))[:],vcat(transpose(y[2*n+1:3*n]),transpose(fill(0,n)))[:],y[3*n+1:3*n+3],gps[1],y[end],gps[3:end]),o,x->sqrt(sum(x)),out,gr_dep)
end

#function to extract maximum liklihood particle
function extract_ml_particle(fit::APMCResult,model)
  temp=fit.populations[model,end]
  md=mode_univar(temp)
  sepmat=broadcast(-,temp,md)
  for i in 1:size(sepmat)[1]
    sepmat[i,:]=sepmat[i,:]./std(sepmat[i,:])
  end
  dists=sqrt.(sum(sepmat.^2,dims=1))
  ind=findmin(dists)
  return((temp[:,ind[2][2]],ind[1]))
end

#error function for poisson simulations: d1 is experimental data, d2 parameters, o the maximum order of moments considered, met the metric for comparison
function rho_general(d1,d2,o=5,met=x->sqrt(sum(x)),out="fit",gr_dep=0)
  nms=size(d1.moms)[1]-1
  #simulate size distribution for 10 cell cycles without gene epxression to approximae steady state before adding gene expression
  lens=nlm(d2[end-4:end],d1.size,10,0.01,true,true)
  #simulate with gene expression for 5 cell cycles
  d2=tvg_lindep(d2[end-4],d2[end-3],d2[end-2],d2[end-1],d2[end],d2[1:2*nms],d2[2*nms+1:4*nms],d2[4*nms+1:6*nms],d2[6*nms+1:8*nms],d1.size,5,0.01,nms,false,lens[:,1],lens[:,2],lens[:,3],lens[:,4],lens[:,5],zeros(nms),ones(nms),zeros(nms),gr_dep)
  #initialise arrays for moments
  d=fill(Inf,size(d1.moms))
  means=Array{Float64}(undef,nms+1)
   vars=Array{Float64}(undef,nms+1)
   # calculate moments from simulations
  for a in 1:nms+1
    means[a]=mean(d2[:,a])
    #vars[a]=mean((d2[:,a]-means[a]).^2)
    vars[a]=1
  end
  a=0
  if nms>1 a=a+1
    for a in 1:nms-1
      d[a,a,a,2,1,1]=means[a]
      if o>1
      for i in 2:o
        d[a,a,a,i+1,1,1]=mean((d2[:,a].-means[a]).^i)/vars[a]^(i/2)
      end
      for b in (a+1):(nms+1)
        for i in 1:(o-1)
          for j in 1:(o-1)
            if (i+j)<(o+1)
              d[a,b,1,i+1,j+1,1]=mean((d2[:,a].-means[a]).^i.*(d2[:,b].-means[b]).^j)/(vars[a]^(i/2)*vars[b]^(j/2))
            end
          end
        end
      end
      end
      if o>2
      for b in (a+1):(nms)
        for c in (b+1):(nms+1)
          for i in 1:(o-2)
            for j in 1:(o-2)
              for k in 1:(o-2)
                if (i+j+k)<(o+1)
                  d[a,b,c,i+1,j+1,k+1]=mean((d2[:,a].-means[a]).^i.*(d2[:,b].-means[b]).^j.*(d2[:,c].-means[c]).^k)/(vars[a]^(i/2)*vars[b]^(j/2)*vars[c]^(k/2))
                end
              end
            end
          end
        end
      end
    end
  end
  end
  a=a+1
  d[a,a,a,2,1,1]=means[a]
  if o>1
  for i in 2:o
    d[a,a,a,i+1,1,1]=mean((d2[:,a].-means[a]).^i)/vars[a]^(i/2)
  end
  for b in (a+1):(nms+1)
    for i in 1:(o-1)
      for j in 1:(o-1)
        if i+j<(o+1)
          d[a,b,1,i+1,j+1,1]=mean((d2[:,a].-means[a]).^i.*(d2[:,b].-means[b]).^j)/(vars[a]^(i/2)*vars[b]^(j/2))
        end
      end
    end
  end
  end
  if (nms>0)
  a=a+1
  d[a,a,a,2,1,1]=means[a]
  if o>1
  for i in 2:o
    d[a,a,a,i+1,1,1]=mean((d2[:,a].-means[a]).^i)/vars[a]^(i/2)
  end
  end
  end
  if out=="moms"
    return(d)
  end
  if out=="fit"
    return(met(((d1.moms[d.!=Inf]-d[d.!=Inf])./(d1.wts[d.!=Inf]).^1).^2)/sqrt(o))
  end
  if out=="raw"
    return(d2)
  end
  println("output not recognised")
end

function rho_cc(d1,d2,o=5,met=x->sqrt(sum(x)),out="fit",gr_dep=0)
  nms=size(d1.moms)[1]-1
  lens=nlm(d2[end-4:end],d1.size,10,0.01,true,true)
  d2=tvg_lindep(d2[end-4],d2[end-3],d2[end-2],d2[end-1],d2[end],d2[1:2*nms],d2[2*nms+1:4*nms],d2[4*nms+1:6*nms],d2[6*nms+1:8*nms],d1.size,5,0.01,nms,false,lens[:,1],lens[:,2],lens[:,3],lens[:,4],lens[:,5],d2[8*nms+1:9*nms],d2[9*nms+1:10*nms],d2[10*nms+1:11*nms],gr_dep)
  d=fill(Inf,size(d1.moms))
  means=Array{Float64}(undef,nms+1)
   vars=Array{Float64}(undef,nms+1)
  for a in 1:nms+1
    means[a]=mean(d2[:,a])
    #vars[a]=mean((d2[:,a]-means[a]).^2)
    vars[a]=1
  end
  a=0
  if nms>1 a=a+1
    for a in 1:nms-1
      d[a,a,a,2,1,1]=means[a]
      if o>1
      for i in 2:o
        d[a,a,a,i+1,1,1]=mean((d2[:,a].-means[a]).^i)/vars[a]^(i/2)
      end
      for b in (a+1):(nms+1)
        for i in 1:(o-1)
          for j in 1:(o-1)
            if (i+j)<(o+1)
              d[a,b,1,i+1,j+1,1]=mean((d2[:,a].-means[a]).^i.*(d2[:,b].-means[b]).^j)/(vars[a]^(i/2)*vars[b]^(j/2))
            end
          end
        end
      end
      end
      if o>2
      for b in (a+1):(nms)
        for c in (b+1):(nms+1)
          for i in 1:(o-2)
            for j in 1:(o-2)
              for k in 1:(o-2)
                if (i+j+k)<(o+1)
                  d[a,b,c,i+1,j+1,k+1]=mean((d2[:,a].-means[a]).^i.*(d2[:,b].-means[b]).^j.*(d2[:,c].-means[c]).^k)/(vars[a]^(i/2)*vars[b]^(j/2)*vars[c]^(k/2))
                end
              end
            end
          end
        end
      end
    end
  end
  end
  a=a+1
  d[a,a,a,2,1,1]=means[a]
  if o>1
  for i in 2:o
    d[a,a,a,i+1,1,1]=mean((d2[:,a].-means[a]).^i)/vars[a]^(i/2)
  end
  for b in (a+1):(nms+1)
    for i in 1:(o-1)
      for j in 1:(o-1)
        if i+j<(o+1)
          d[a,b,1,i+1,j+1,1]=mean((d2[:,a].-means[a]).^i.*(d2[:,b].-means[b]).^j)/(vars[a]^(i/2)*vars[b]^(j/2))
        end
      end
    end
  end
  end
  a=a+1
  d[a,a,a,2,1,1]=means[a]
  if o>1
  for i in 2:o
    d[a,a,a,i+1,1,1]=mean((d2[:,a].-means[a]).^i)/vars[a]^(i/2)
  end
  end
  if out=="moms"
    return(d)
  end
  if out=="fit"
    return(met(((d1.moms[d.!=Inf]-d[d.!=Inf])./(d1.wts[d.!=Inf]).^1).^2)/sqrt(o))
  end
  if out=="raw"
    return(d2)
  end
  println("output not recognised")
end

function rho_bursty(d1,d2,o=5,met=x->sqrt(sum(x)),out="fit",gr_dep=0)
  nms=size(d1.moms)[1]-1
  lens=nlm(d2[end-4:end],d1.size,10,0.01,true,true)
  d2=bursty_lindep(d2[end-4],d2[end-3],d2[end-2],d2[end-1],d2[end],d2[1:2*nms],d2[2*nms+1:4*nms],d2[4*nms+1:6*nms],d1.size,5,0.01,nms,false,lens[:,1],lens[:,2],lens[:,3],lens[:,4],lens[:,5],zeros(nms),ones(nms),zeros(nms),gr_dep)
  d=fill(Inf,size(d1.moms))
  means=Array{Float64}(undef,nms+1)
   vars=Array{Float64}(undef,nms+1)
  for a in 1:nms+1
    means[a]=mean(d2[:,a])
    #vars[a]=mean((d2[:,a]-means[a]).^2)
    vars[a]=1
  end
  a=0
  if nms>1 a=a+1
    for a in 1:nms-1
      d[a,a,a,2,1,1]=means[a]
      if o>1
      for i in 2:o
        d[a,a,a,i+1,1,1]=mean((d2[:,a].-means[a]).^i)/vars[a]^(i/2)
      end
      for b in (a+1):(nms+1)
        for i in 1:(o-1)
          for j in 1:(o-1)
            if (i+j)<(o+1)
              d[a,b,1,i+1,j+1,1]=mean((d2[:,a].-means[a]).^i.*(d2[:,b].-means[b]).^j)/(vars[a]^(i/2)*vars[b]^(j/2))
            end
          end
        end
      end
      end
      if o>2
      for b in (a+1):(nms)
        for c in (b+1):(nms+1)
          for i in 1:(o-2)
            for j in 1:(o-2)
              for k in 1:(o-2)
                if (i+j+k)<(o+1)
                  d[a,b,c,i+1,j+1,k+1]=mean((d2[:,a].-means[a]).^i.*(d2[:,b].-means[b]).^j.*(d2[:,c].-means[c]).^k)/(vars[a]^(i/2)*vars[b]^(j/2)*vars[c]^(k/2))
                end
              end
            end
          end
        end
      end
    end
  end
  end
  a=a+1
  d[a,a,a,2,1,1]=means[a]
  if o>1
  for i in 2:o
    d[a,a,a,i+1,1,1]=mean((d2[:,a].-means[a]).^i)/vars[a]^(i/2)
  end
  for b in (a+1):(nms+1)
    for i in 1:(o-1)
      for j in 1:(o-1)
        if i+j<(o+1)
          d[a,b,1,i+1,j+1,1]=mean((d2[:,a].-means[a]).^i.*(d2[:,b].-means[b]).^j)/(vars[a]^(i/2)*vars[b]^(j/2))
        end
      end
    end
  end
  end
  a=a+1
  d[a,a,a,2,1,1]=means[a]
  if o>1
  for i in 2:o
    d[a,a,a,i+1,1,1]=mean((d2[:,a].-means[a]).^i)/vars[a]^(i/2)
  end
  end
  if out=="moms"
    return(d)
  end
  if out=="fit"
    return(met(((d1.moms[d.!=Inf]-d[d.!=Inf])./(d1.wts[d.!=Inf]).^1).^2)/sqrt(o))
  end
  if out=="raw"
    return(d2)
  end
  println("output not recognised")
end

function rho_bursty_cc(d1,d2,o=5,met=x->sqrt(sum(x)),out="fit",gr_dep=0)
  nms=size(d1.moms)[1]-1
  lens=nlm(d2[end-4:end],d1.size,10,0.01,true,true)
  d2=bursty_lindep(d2[end-4],d2[end-3],d2[end-2],d2[end-1],d2[end],d2[1:2*nms],d2[2*nms+1:4*nms],d2[4*nms+1:6*nms],d1.size,5,0.01,nms,false,lens[:,1],lens[:,2],lens[:,3],lens[:,4],lens[:,5],d2[6*nms+1:8*nms],d2[8*nms+1:10*nms],d2[10*nms+1:12*nms],gr_dep)
  d=fill(Inf,size(d1.moms))
  means=Array{Float64}(undef,nms+1)
   vars=Array{Float64}(undef,nms+1)
  for a in 1:nms+1
    means[a]=mean(d2[:,a])
    #vars[a]=mean((d2[:,a]-means[a]).^2)
    vars[a]=1
  end
  a=0
  if nms>1 a=a+1
    for a in 1:nms-1
      d[a,a,a,2,1,1]=means[a]
      if o>1
      for i in 2:o
        d[a,a,a,i+1,1,1]=mean((d2[:,a].-means[a]).^i)/vars[a]^(i/2)
      end
      for b in (a+1):(nms+1)
        for i in 1:(o-1)
          for j in 1:(o-1)
            if (i+j)<(o+1)
              d[a,b,1,i+1,j+1,1]=mean((d2[:,a].-means[a]).^i.*(d2[:,b].-means[b]).^j)/(vars[a]^(i/2)*vars[b]^(j/2))
            end
          end
        end
      end
      end
      if o>2
      for b in (a+1):(nms)
        for c in (b+1):(nms+1)
          for i in 1:(o-2)
            for j in 1:(o-2)
              for k in 1:(o-2)
                if (i+j+k)<(o+1)
                  d[a,b,c,i+1,j+1,k+1]=mean((d2[:,a].-means[a]).^i.*(d2[:,b].-means[b]).^j.*(d2[:,c].-means[c]).^k)/(vars[a]^(i/2)*vars[b]^(j/2)*vars[c]^(k/2))
                end
              end
            end
          end
        end
      end
    end
  end
  end
  a=a+1
  d[a,a,a,2,1,1]=means[a]
  if o>1
  for i in 2:o
    d[a,a,a,i+1,1,1]=mean((d2[:,a].-means[a]).^i)/vars[a]^(i/2)
  end
  for b in (a+1):(nms+1)
    for i in 1:(o-1)
      for j in 1:(o-1)
        if i+j<(o+1)
          d[a,b,1,i+1,j+1,1]=mean((d2[:,a].-means[a]).^i.*(d2[:,b].-means[b]).^j)/(vars[a]^(i/2)*vars[b]^(j/2))
        end
      end
    end
  end
  end
  a=a+1
  d[a,a,a,2,1,1]=means[a]
  if o>1
  for i in 2:o
    d[a,a,a,i+1,1,1]=mean((d2[:,a].-means[a]).^i)/vars[a]^(i/2)
  end
  end
  if out=="moms"
    return(d)
  end
  if out=="fit"
    return(met(((d1.moms[d.!=Inf]-d[d.!=Inf])./(d1.wts[d.!=Inf]).^1).^2)/sqrt(o))
  end
  if out=="raw"
    return(d2)
  end
  println("output not recognised")
end

#Noisy Linear Map code
function nlm_init(params,N)
    a=params[1]
    b=params[2]
    eps1=params[3]
    eps2=params[4]
    eps3=params[5]
    raw=Array{Float64}(undef,N)
    d1=Normal(0.0,eps1^2)
    d2=Normal(0.5,eps2)
    if eps3!=0
      d3=TruncatedNormal(1.0,eps3,0.0,Inf)
    end
    d4=Normal(0.0,eps1/(1-(a^2)/4))
    starts=max.(rand(d2,N).*(fill(2*b/(2-a),N).+rand(d4,N)),fill(0.25*b/(2.0-a),N))
    fins=max.(a*starts.+b.+rand(d1,N),starts)
    t0=2 .^rand(Uniform(0,1),N).-fill(1.0,N)
    raw[:]=starts.+(fins.-starts).*t0
    if eps3==0
      grs=fill(1,N)
    else
      grs=rand(d3,N)
    end
    div_waits=log.(2,fins./starts).*(1 .-t0)./grs
    return starts, fins, raw, grs, div_waits
end

function nlm_div!(i,t,dt,raw,div_waits,starts,fins,grs,gf,params,N,unbiased)
  a=params[1]
  b=params[2]
  eps1=params[3]
  eps2=params[4]
  eps3=params[5]
  d1=Normal(0.0,eps1^2)
  d2=Normal(0.5,eps2)
  if eps3!=0
    d3=TruncatedNormal(1.0,eps3,0.0,Inf)
  end
  div=min(max(rand(d2),0),1)
  if unbiased
    dtr=sample(1:(N+1))
  else
    dtr=i
  end
  if dtr<(N+1)
    starts[dtr]=fins[i]*(1-div)
    if eps3!=0
      grs[dtr]=rand(d3)
      gf[dtr]=2^(grs[dtr]*dt)
    end
    fins[dtr]=max(a*starts[dtr]+b+rand(d1),starts[dtr])
    raw[dtr]=starts[dtr]*2^(-grs[dtr]*div_waits[i])
    div_waits[dtr]=div_waits[i]+log(2,fins[dtr]/starts[dtr])/grs[dtr]
  end
  if dtr!=i
    starts[i]=fins[i]*div
    if eps3!=0
      grs[i]=rand(d3)
      gf[i]=2^(grs[i]*dt)
    end
    raw[i]=starts[i]*2^(-grs[i]*div_waits[i])
    fins[i]=max(a*starts[i]+b+rand(d1),starts[i])
    div_waits[i]=div_waits[i]+log(2,fins[i]/starts[i])/grs[i]
  end
  return(div)
end

function nlm(params,N,tf,dt,unbiased=true,times=false)
  starts, fins, raw, grs, div_waits =nlm_init(params,N)
  gf=2 .^(grs.*dt)
  t=0
  while t<tf
    t=t+dt
    for i in 1:N
      div_waits[i]=div_waits[i]-dt
      if div_waits[i]<=0
        nlm_div!(i,t,dt,raw,div_waits,starts,fins,grs,gf,params,N,unbiased)
      else
        raw[i]=min(raw[i]*gf[i],fins[i])
      end
    end
  end
  if times
    return hcat(starts,fins,div_waits,grs,raw)
  else
    return raw
  end
end

# time varying Gillespie code for non bursty models
function tvg_lindep(a,b,eps1,eps2,eps3,kon,koff,v,d,N,tf,dt,nms=div(length(v),2),scratch=true,starts=zeros(N),fins=zeros(N),div_waits=zeros(N),grs=zeros(N),vals=zeros(N),ton=zeros(nms),toff=ones(nms),resid=zeros(nms),gr_dep=0)
  raw=Array{Float64}(undef,N,nms+1)
  parts=Array{Int}(undef,nms)
  d1=Normal(0.0,eps1)
  d2=TruncatedNormal(0.5,eps2,0,1)
  d3=Uniform()
  d4=Exponential()
  as=Array{Float64}(undef,4*nms)
  if eps3!=0
  d5=TruncatedNormal(1.0,eps3,0.0,Inf)
  end
  if scratch
    if eps3!=0
        grs=rand(d5,N)
    else
      grs=ones(N)
    end
  init=Uniform(0.5*b/(2-a),2*b/(2-a))
  starts=rand(init,N)
  fins=max.(a.*starts.+b.+rand(d1,N),starts)
  t0=rand(d3,N)
  raw[:,nms+1]=starts.+(fins.-starts).*t0
  div_waits=log.(2,fins./starts).*(1 .-t0)./grs
  prog= log.(2,fins./starts).*(1 .-t0)./grs
  else
  raw[:,nms+1]=vals
    prog=fill(Inf,N)
  end
  gf=2  .^(dt.*grs)
  switch=ones(N,nms)
  raw[:,1:nms]=zeros(N,nms)
  t=0
  while t<tf
    for i in 1:N
      k=min(0.5*(1+gf[i])*raw[i,nms+1],100)
      if gr_dep==0
          k2=1
      elseif gr_dep==1
          k2=grs[i]
      else
          println("invalid gr_dep")
          break
      end
      for j in 1:nms
          if j==1
      as[4*(j-1)+1]=max((v[2*(j-1)+1]+v[2*(j-1)+2]*k)*(switch[i,j])*k2*((((1-div_waits[i]/prog[i])>=ton[j])&((1-div_waits[i]/prog[i])<toff[j]))+resid[j]*(((1-div_waits[i]/prog[i])<ton[j])|((1-div_waits[i]/prog[i])>=toff[j]))),0.0)
          else
      as[4*(j-1)+1]=as[4*(j-1)]+max((v[2*(j-1)+1]+v[2*(j-1)+2]*k)*(switch[i,j])*k2*((((1-div_waits[i]/prog[i])>=ton[j])&((1-div_waits[i]/prog[i])<toff[j]))+resid[j]*(((1-div_waits[i]/prog[i])<ton[j])|((1-div_waits[i]/prog[i])>=toff[j]))),0.0)
          end
      as[4*(j-1)+2]=as[4*(j-1)+1]+abs((d[2*(j-1)+1]+d[2*(j-1)+2]/k)*raw[i,j])
      as[4*(j-1)+3]=as[4*(j-1)+2]+max((kon[2*(j-1)+1]+kon[2*(j-1)+2]*k)*(switch[i,j]==resid),0.0)
      as[4*(j-1)+4]=as[4*(j-1)+3]+max((koff[2*(j-1)+1]+koff[2*(j-1)+2]*k)*(switch[i,j]==1),0.0)
        end
      a0=as[end]
      tg=rand(d4)/a0
      tgs=tg
      while (tgs<dt)&(tgs<(div_waits[i]))
        u=a0*rand(d3)
        mrna=1
        while u>as[(4*mrna)]
          mrna+=1
        end
        if u<=as[4*(mrna-1)+1]
          raw[i,mrna]+=1
        else
        if u<=as[4*(mrna-1)+2]
            raw[i,mrna]+=-1
          else
        if u<=as[4*(mrna-1)+3]
          switch[i,mrna]=1
        else
          switch[i,mrna]=resid[mrna]
        end
          end
        end
        for j in 1:nms
            if j==1
        as[4*(j-1)+1]=max((v[2*(j-1)+1]+v[2*(j-1)+2]*k)*(switch[i,j])*k2*((((1-div_waits[i]/prog[i])>=ton[j])&((1-div_waits[i]/prog[i])<toff[j]))+resid[j]*(((1-div_waits[i]/prog[i])<ton[j])|((1-div_waits[i]/prog[i])>=toff[j]))),0.0)
            else
        as[4*(j-1)+1]=as[4*(j-1)]+max((v[2*(j-1)+1]+v[2*(j-1)+2]*k)*(switch[i,j])*k2*((((1-div_waits[i]/prog[i])>=ton[j])&((1-div_waits[i]/prog[i])<toff[j]))+resid[j]*(((1-div_waits[i]/prog[i])<ton[j])|((1-div_waits[i]/prog[i])>=toff[j]))),0.0)
            end
        as[4*(j-1)+2]=as[4*(j-1)+1]+abs((d[2*(j-1)+1]+d[2*(j-1)+2]/k)*raw[i,j])
        as[4*(j-1)+3]=as[4*(j-1)+2]+max((kon[2*(j-1)+1]+kon[2*(j-1)+2]*k)*(switch[i,j]==resid),0.0)
        as[4*(j-1)+4]=as[4*(j-1)+3]+max((koff[2*(j-1)+1]+koff[2*(j-1)+2]*k)*(switch[i,j]==1),0.0)
          end
        a0=as[end]
        tg=rand(d4)/a0
        tgs+=tg
      end
      div_waits[i]=div_waits[i]-dt
      if(div_waits[i]<0)
        div=rand(d2)
        for m in 1:nms
        parts[m]=rand(Binomial(round(Int,raw[i,m]),div))
        end
        dtr=sample(1:(N+1))
        if dtr<(N+1)
       raw[dtr,1:nms]=raw[i,1:nms]-parts
        switch[dtr,1:nms]=switch[i,1:nms]
        starts[dtr]=fins[i]*(1-div)
          if eps3!=0
        grs[dtr]=rand(d5)
            gf[dtr]=2^(grs[dtr]*dt)
          end
        fins[dtr]=max(a*starts[dtr]+b+rand(d1),starts[dtr])
        raw[dtr,nms+1]=starts[dtr]*2^(-grs[dtr]*div_waits[i])
        div_waits[dtr]=div_waits[i]+log(2,fins[dtr]/starts[dtr])/grs[dtr]
          prog[dtr]=log(2,fins[dtr]/starts[dtr])/grs[dtr]
        end
        if dtr!=i
        starts[i]=fins[i]*div
          if eps3!=0
        grs[i]=rand(d5)
            gf[i]=2^(grs[i]*dt)
          end
        raw[i,nms+1]=starts[i]*2^(-grs[i]*div_waits[i])
        raw[i,1:nms]=transpose(parts)
        fins[i]=max(a*starts[i]+b+rand(d1),starts[i])
        div_waits[i]=div_waits[i]+log(2,fins[i]/starts[i])/grs[i]
          prog[i]=log(2,fins[i]/starts[i])/grs[i]
        end
      else
        raw[i,nms+1]=raw[i,nms+1]*gf[i]
      end
  end
       t=t+dt
  end
  return(raw)
end

#time varying Gillespie code for bursty models (geometric distribution for burst sizes)
function bursty_lindep(a,b,eps1,eps2,eps3,kon,bs,d,N,tf,dt,nms=1,scratch=true,starts=zeros(N),fins=zeros(N),div_waits=zeros(N),grs=zeros(N),vals=zeros(N),ton=zeros(nms),toff=ones(nms),resid=zeros(nms),gr_dep=0)
  raw=Array{Float64}(undef,N,nms+1)
    parts=Array{Int}(undef,nms)
  d1=Normal(0.0,eps1)
  d2=TruncatedNormal(0.5,eps2,0,1)
  d3=Uniform()
  d4=Exponential()
  as=Array{Float64}(undef,2*nms)
  if eps3!=0
  d5=TruncatedNormal(1.0,eps3,0.0,Inf)
  end
  if scratch
    if eps3!=0
        grs=rand(d5,N)
    else
      grs=ones(N)
    end
  init=Uniform(0.5*b/(2-a),2*b/(2-a))
  starts=rand(init,N)
  fins=max.(a*starts.+b.+rand(d1,N),starts)
  t0=rand(d3,N)
  raw[:,nms+1]=starts.+(fins.-starts).*t0
  div_waits=log.(2,fins./starts).*(1 .-t0)./grs
  prog= log.(2,fins./starts).*(1 .-t0)./grs
  else
  raw[:,nms+1]=vals
    prog=fill(Inf,N)
  end
  gf=2 .^(dt.*grs)
  raw[:,1:nms]=zeros(N,nms)
  t=0.0
  while t<tf
    for i in 1:N
      k=0.5*(1+gf[i])*raw[i,nms+1]
      k2=1+(grs[i]-1)*gr_dep
      for j in 1:nms
      as[2*(j-1)+1]=min(max((kon[2*(j-1)+1]+kon[2*(j-1)+2]*k)*((((1-div_waits[i]/prog[i])>=ton[j])&((1-div_waits[i]/prog[i])<toff[j]))+resid[j]*(((1-div_waits[i]/prog[i])<ton[j])|((1-div_waits[i]/prog[i])>=toff[j]))),0.0),1000)
      as[2*(j-1)+2]=min(max((d[2*(j-1)+1]+d[2*(j-1)+2]/k)*raw[i,j],0.0),1000)
      end
      a0=sum(as)
      tg=rand(d4)/a0
      tgs=tg
      while (tgs<dt)&(tgs<(div_waits[i]-dt))
        u=a0*rand(d3)
        mrna=1
        s=0
        s2=s+sum(as[(2*(mrna-1)+1):(2*mrna)])
        while u>s2
          mrna+=1
          s=s2
          s2+=sum(as[(2*(mrna-1)+1):(2*mrna)])
        end
        s+=as[2*(mrna-1)+1]
        if u<=s
          #println((fld(log(rand(d3)),log(1-min(max(k2*(bs[2*(mrna-1)+1]+bs[2*(mrna-1)+2]*k)/(1+k2*bs[2*(mrna-1)+1]+bs[2*(mrna-1)+2]*k),0.0),1.0)))+1))
          #println(fld(log(rand(d3)),log(1-min(max(1/(k2*(bs[2*(mrna-1)+1]+bs[2*(mrna-1)+2]*k)),0.0),1.0)))+1)
         raw[i,mrna]+=(fld(log(rand(d3)),log(1-min(max(1/(k2*(bs[2*(mrna-1)+1]+bs[2*(mrna-1)+2]*k)+1),0.001),1)))+0)
        else
         raw[i,mrna]+=-1
        end
      for j in 1:nms
          as[2*(j-1)+1]=min(max((kon[2*(j-1)+1]+kon[2*(j-1)+2]*k)*((((1-div_waits[i]/prog[i])>=ton[j])&((1-div_waits[i]/prog[i])<toff[j]))+resid[j]*(((1-div_waits[i]/prog[i])<ton[j])|((1-div_waits[i]/prog[i])>=toff[j]))),0.0),1000)
      as[2*(j-1)+2]=min(max((d[2*(j-1)+1]+d[2*(j-1)+2]/k)*raw[i,j],0.0),1000)
      end
      a0=sum(as)
        tg=rand(d4)/a0
        tgs+=tg
      end
      div_waits[i]=div_waits[i]-dt
      if(div_waits[i]<0)
        div=rand(d2)
        for m in 1:nms
        parts[m]=rand(Binomial(round(Int,raw[i,m]),div))
        end
        dtr=sample(1:(N+1))
        if dtr<(N+1)
       raw[dtr,1:nms]=raw[i,1:nms]-transpose(parts)
        starts[dtr]=fins[i]*(1-div)
          if eps3!=0
        grs[dtr]=rand(d5)
            gf[dtr]=2^(grs[dtr]*dt)
          end
        fins[dtr]=max(a*starts[dtr]+b+rand(d1),starts[dtr])
        raw[dtr,nms+1]=starts[dtr]*2^(-grs[dtr]*div_waits[i])
        div_waits[dtr]=div_waits[i]+log(2,fins[dtr]/starts[dtr])/grs[dtr]
        prog[dtr]=log(2,fins[dtr]/starts[dtr])/grs[dtr]
      end
      if dtr!=i
      starts[i]=fins[i]*div
        if eps3!=0
      grs[i]=rand(d5)
          gf[i]=2^(grs[i]*dt)
        end
      raw[i,nms+1]=starts[i]*2^(-grs[i]*div_waits[i])
      raw[i,1:nms]=transpose(parts)
      fins[i]=max(a*starts[i]+b+rand(d1),starts[i])
      div_waits[i]=div_waits[i]+log(2,fins[i]/starts[i])/grs[i]
        prog[i]=log(2,fins[i]/starts[i])/grs[i]
      end
      else
        raw[i,nms+1]=raw[i,nms+1]*gf[i]
      end
  end
       t=t+dt
  end
  return(raw)
end

#set up parallel backend
addprocs(3)
@everywhere include(joinpath(homedir(),"Box","Dropbox","FishProjectR","julia","juliarc.jl"))

#size distribution prefits
cd("C:\\Users\\ab5910\\Box\\Dropbox\\FishProjectR\\Julia")
ref=readdir("lengths")
lens=Array{DataFrame}(undef,length(ref))
for i in 1:length(ref)
    lens[i]=CSV.read(joinpath(homedir(),"Box Sync","Dropbox","FishProjectR","Julia","lengths",ref[i]),header=1)
end
lens_fits=Array{ABCfit}(undef,length(ref))
for i in 1:length(ref)
    println(i)
    lens_fits[i]=APMC_KDE(2000,lens[i][:,end]+rand(Normal(0,0.00001),length(lens[i][:,end])),Vector[model_lens,model_lens_SE],[(x,y)->rho_lens(x,y),(x,y)->rho_lens_SE(x,y)])
end
nlm_bases=map(extract_ml_particle,lens_fits,fill(1,11))

#fig2 fits

ref=vcat(readdir("FIT_DATA"),readdir("SIB1_DATA"))
expds_fig2=Array{gsummary}(undef,length(ref))
for i in 1:41
    println(i)
    expds_fig2[i]=calc_moms(CSV.read(joinpath(homedir(),"Box Sync","Dropbox","FishProjectR","Julia","FIT_DATA",ref[i]),header=1)[:,3:4])
end
for i in 42:50
    println(i)
    expds_fig2[i]=calc_moms(CSV.read(joinpath(homedir(),"Box Sync","Dropbox","FishProjectR","Julia","SIB1_DATA",ref[i]),header=1)[:,3:4])
end

#wild type YE constituitive genes
ye_wt_const=[7,13,18,21,24,27,30]
mods=[model_pois_preset(1),model_pois_v_preset(1),model_bursty_preset(1),model_bursty_v_preset(1),model_bursty_kon_preset(1)]
rhos=[f_pois_preset(1,nlm_bases[4][1],0),f_pois_v_preset(1,nlm_bases[4][1],0),f_bursty_preset(1,nlm_bases[4][1],0),f_bursty_v_preset(1,nlm_bases[4][1],0),f_bursty_kon_preset(1,nlm_bases[4][1],0)]
for i in ye_wt_const
    println(ref[i])
    fits_fig2[i]=APMC_KDE(10000,expds_fig2[i],mods,rhos,paccmin=0.05)
end

#wild type EMM constituitive genes
emm_wt=occursin.("EMM",ref).&occursin.("972h-",ref)
emm_wt=findall(emm_wt)
mods=[model_pois_preset(1),model_pois_v_preset(1),model_bursty_preset(1),model_bursty_v_preset(1),model_bursty_kon_preset(1)]
rhos=[f_pois_preset(1,nlm_bases[2][1],0),f_pois_v_preset(1,nlm_bases[2][1],0),f_bursty_preset(1,nlm_bases[2][1],0),f_bursty_v_preset(1,nlm_bases[2][1],0),f_bursty_kon_preset(1,nlm_bases[2][1],0)]
for i in emm_wt
    println(ref[i])
    fits_fig2[i]=APMC_KDE(10000,expds_fig2[i],mods,rhos,paccmin=0.05)
end

#pom1
misc=.!(occursin.("972h-",ref).|occursin.("wee1",ref).|occursin.("cdc25",ref))
misc=findall(misc)
pom1=misc[2]
mods=[model_pois_preset(1),model_pois_v_preset(1),model_bursty_preset(1),model_bursty_v_preset(1),model_bursty_kon_preset(1)]
rhos=[f_pois_preset(1,nlm_bases[8][1],0),f_pois_v_preset(1,nlm_bases[8][1],0),f_bursty_preset(1,nlm_bases[8][1],0),f_bursty_v_preset(1,nlm_bases[8][1],0),f_bursty_kon_preset(1,nlm_bases[8][1],0)]
fits_fig2[pom1]=APMC_KDE(10000,expds_fig2[pom1],mods,rhos,paccmin=0.05)

#wild type YE cell cycle genes
ye_wt_cc=[1,4,10]
mods=[model_pois_cc_preset(1),model_pois_cc_v_preset(1),model_bursty_cc_preset(1),model_bursty_cc_v_preset(1),model_bursty_cc_kon_preset(1)]
rhos=[f_pois_cc_preset(1,nlm_bases[4][1],0),f_pois_cc_v_preset(1,nlm_bases[4][1],0),f_bursty_cc_preset(1,nlm_bases[4][1],0),f_bursty_cc_v_preset(1,nlm_bases[4][1],0),f_bursty_cc_kon_preset(1,nlm_bases[4][1],0)]
for i in ye_wt_cc
    println(ref[i])
    fits_fig2[i]=APMC_KDE(10000,expds_fig2[i],mods,rhos,paccmin=0.05)
end

#sib1
ye_wt_sib1=[42,43,44]
mods=[model_pois_preset(1),model_pois_v_preset(1),model_bursty_preset(1),model_bursty_v_preset(1),model_bursty_kon_preset(1)]
rhos=[f_pois_preset(1,nlm_bases[4][1],0),f_pois_v_preset(1,nlm_bases[4][1],0),f_bursty_preset(1,nlm_bases[4][1],0),f_bursty_v_preset(1,nlm_bases[4][1],0),f_bursty_kon_preset(1,nlm_bases[4][1],0)]
for i in ye_wt_sib1
    println(ref[i])
    fits_fig2[i]=APMC_KDE(10000,expds_fig2[i],mods,rhos,paccmin=0.05)
end
