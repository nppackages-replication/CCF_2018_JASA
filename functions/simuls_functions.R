mu.fun <- function(m,kern) integrate(function(x) x^(m)*kern(x),-Inf,Inf)$value
vartheta.fun <- function(p,kern) integrate(function(x) (kern(x))^p,-Inf,Inf)$value

K.0.2.mv = function(x) {(abs(x)<=1)*(1/2)}
K.0.4.mv = function(x) {(abs(x)<=1)*(3/8)*(3-5*x^2)}
K.2.4.mv = function(x) {(abs(x)<=1)*(15/4)*(3*x^2-1)}

K.0.2.ms = function(x) {(abs(x)<=1)*(3/4)*(1-x^2)}
K.0.4.ms = function(x) {(abs(x)<=1)*(15/32)*(7*x^4-10*x^2+3)}
K.2.4.ms = function(x) {(abs(x)<=1)*(105/16)*(6*x^2-5*x^4-1)}

#M.epan  = function(u,rho,mu) {k.epan(u)-rho^(1+2)*L.2(rho*u)*mu}

M.unif.mv = function(x) { K.0.2.mv(x) - 0.5*rho^(1+r)*K.2.4.mv(rho*x)*mu.fun(2,K.0.2.mv) }
M.epan.ms = function(x) { K.0.2.ms(x) - 0.5*rho^(1+r)*K.2.4.ms(rho*x)*mu.fun(2,K.0.2.ms) }

kd.K.fun = function(x,v,r,kernel){
  if (v==2) {
    if (kernel=="gau"){
      if (r==0) k = function(u)            dnorm(u)
      if (r==2) k = function(u)     (u^2-1)*dnorm(u)
      if (r==4) k = function(u) (u^4-6*u^2+3)*dnorm(u)
    }
    if (kernel=="uni"){
      if (r==0) k = function(u)  0.5*(abs(u)<=1)
    }
    if (kernel=="epa"){
      if (r==0) k = function(u)  0.75*(1-u^2)*(abs(u)<=1)
    }
  }
  if (v==4) {
    if (kernel=="uni"){
      if (r==0) k = function(u)    (abs(u)<=1)* 3*(-5*u^2+3)/8
      if (r==2) k = function(u)    (abs(u)<=1)*15*(3*u^2-1)/4
    }
    if (kernel=="epa"){
      if (r==0) k = function(u)   (abs(u)<=1)*(15/32)*(7*u^4-10*u^2+3)
      if (r==2) k = function(u)   (abs(u)<=1)*(105/16)*(6*u^2-5*u^4-1)
    }
  }
  if (v==6) {
    if (kernel=="uni"){
      if (r==0) k = function(u)    (abs(u)<=1)*15*(63*u^4-70*u^2+15)/128
      if (r==2) k = function(u)    (abs(u)<=1)*105*(-45*u^4+42*u^2-5)/32
    }
    if (kernel=="epa"){
      if (r==0) k = function(u)   (abs(u)<=1)*(35/256)*(-99*u^6+189*x^4-105*x^2+15)
      if (r==2) k = function(u)   (abs(u)<=1)*(315/64)*(77*x^6-135*x^4+63*x^2-5)
    }
  }
  Kx = k(x)
  k.v = integrate(function(u) (-1)^v*u^(v)*k(u)/factorial(v), -Inf,Inf)$value
  R.v = integrate(function(u) (k(u))^2,   -Inf,Inf)$value
  out=list(Kx=Kx, k.v=k.v, R.v=R.v)
  return(out)
}


if (model==1) g.0.fun = function(x) {sin(4*x) + 2*exp(-4*16*x^2)}
if (model==2) g.0.fun = function(x) {2*x+2*exp(-4*16*x^2)} 
if (model==3) g.0.fun = function(x) {0.3*exp(-4*(2*x+1)^2) + 0.7*exp(-16*(2*x-1)^2) }
if (model==4) g.0.fun = function(x) {x+5*dnorm(10*x)}
if (model==5) g.0.fun = function(x) {sin(3*pi*x/2)/(1+18*x^2*(sign(x)+1))}
if (model==6) g.0.fun = function(x) {sin(  pi*x/2)/(1+ 2*x^2*(sign(x)+1))}

x.gen   = function(n) {runif(n, min=-1, max=1)}  
fx      = function(x) {dunif(x, min=-1, max=1)}
u.rnd   = function(n) {rnorm(n,0,s2.pob)}

m1 = function(i,j,k)   integrate(function(x) x^i*x^j*k(x),0,Inf)$value
m2 = function(i,j,k)   integrate(function(x) x^i*x^j*(k(x))^2,0,Inf)$value
m3 = function(i,j,k,a) integrate(function(x) x^i*(a*x)^j*k(x)*k(a*x),0,Inf)$value
GAMMA = function(p,k) {out=matrix(NA,p+1,p+1); for (i in 0:p) for (j in 0:p) out[i+1,j+1]=m1(i,j,k); return(out)}
NU    = function(p,k) {out=matrix(NA,p+1,1); for (i in 0:p) out[i+1,1]=m1(i,p+1,k); return(out)}
PSI   = function(p,k) {out=matrix(NA,p+1,p+1); for (i in 0:p) for (j in 0:p) out[i+1,j+1]=m2(i,j,k); return(out)}
B.lp  = function(p,k) {out=solve(GAMMA(p,k))%*%NU(p,k); out[1]}
C1.fun = function(p0,v,K) {
  S.inv = solve(GAMMA(p0,K)) 
  C1 = (S.inv%*%NU(p0,K))[v+1]
  return(C1)
}
C2.fun = function(p0,v,K) {
  S.inv = solve(GAMMA(p0,K)) 
  C2 = (S.inv%*%PSI(p0,K)%*%S.inv)[v+1,v+1]
  return(C2)
}

h.fun = function(V,B,p,v,N) {(V/(N*B))^(1/(2*p+3))}
N.fun = function(v,V)       {(2*v+1)*V}
D.fun = function(p,v,B,K)   {2*(p+1-v)*(C1.fun(p,v,K)*B)^2}

qrreg = function(x,y,w,s2=0,var.comp=TRUE, ...) {
  M.X = sqrt(w)*x
  X.M.X_inv = qrXXinv(M.X) 
  X.M.Y = crossprod(M.X,sqrt(w)*y)
  beta.hat = X.M.X_inv%*%X.M.Y
  Psi.hat=Sigma.hat=0
  if (var.comp==TRUE) {
    Psi.hat = crossprod((w*s2*w)*x,x)
    Sigma.hat = crossprod(Psi.hat%*%X.M.X_inv,X.M.X_inv)
  }
  output = list(X.M.X_inv=X.M.X_inv, X.M.Y=X.M.Y, beta.hat=beta.hat, Psi.hat=Psi.hat, Sigma.hat=Sigma.hat)
  return(output)
}


m.hat = function(y, x, c, p, h, kernel, deriv) {
  n=length(x)
  w = W.fun((x-c)/h, kernel)/h
  ind = w>0
  eN = sum(ind)
  eY = y[ind]
  eX = x[ind]
  eW = w[ind]
  R = matrix(NA,eN,(p+1))
  for (j in 1:(p+1))  R[,j] = (eX-c)^(j-1)
  invG  = qrXXinv((sqrt(eW)*R))
  tau   = factorial(deriv)*(invG%*%crossprod(R*eW,eY))[(deriv+1),1]
  return(tau)
}

hh.fun = function(y,x,c,h) {
  B = 500
  n=length(y)
  kappa = 2/3
  alpha.0 = 0.05
  z.alpha2 = qnorm(1-alpha.0/2)
  cov = 1-alpha.0
  alpha.seq = seq(0.01,0.09, 0.002)
  x.seq = seq(-0.9,0.9,0.1)
  vareps = 0.1
  nx.seq=length(x.seq)
  gx.xseq   = f.xseq = se.xseq = quant.xseq = gx.hh.b = se.hh.b = 0
  gx.xseq.b = se.xseq.b = matrix(NA,B,nx.seq)
  
  lfit = locfit(y~lp(x, deg = 1))
  y.hat = fitted(lfit)
  res = y - y.hat
  res.hat = res - mean(res)
  s2.hat = mean(res.hat^2)

  gx.hh = m.hat(y=y, x=x, c=c, p=p, h=h, kernel=kernel, deriv=0) 
  h1 = bw.SJ(x)
  f.hh = mean(W.fun((x-c)/h,kernel))/h1
  se.hh = sqrt((kappa*s2.hat)/(n*h*f.hh))
  
  for (x0 in 1:nx.seq) {
    gx.xseq[x0] = m.hat(y=y, x=x, c=x.seq[x0], p=p, h=h, kernel=kernel, deriv=0) 
    f.xseq[x0]   = mean(W.fun((x-x.seq[x0])/h,kernel))/h1
    se.xseq[x0]  = sqrt((kappa*s2.hat)/(n*h*f.xseq[x0]))
  }
  
  #### Bootstrap
  for (b in 1:B) {
    b_ind = sample(n,replace = TRUE)
    res.b = res.hat[b_ind]
    y.b = y.hat + res.b
    lfit.b = locfit(y.b~lp(x, deg = 1))
    y.hat.b = fitted(lfit.b)
    res.b = y.b - y.hat.b
    res.hat.b = res.b - mean(res.b)
    s2.hat.b = mean(res.hat.b^2)
  
    for (x0 in 1:nx.seq) {
      gx.xseq.b[b,x0] = m.hat(y=y.b, x=x, c=x.seq[x0], p=p, h=h, kernel=kernel, deriv=0) 
      se.xseq.b[b,x0] = sqrt((kappa*s2.hat.b)/(n*h*f.xseq[x0]))
    }
    
  }
  
  alpha.hat = 0
  for (x0 in 1:nx.seq) {
    pi.hat = 0
    for (k in 1:length(alpha.seq)) {
      quant.k = qnorm(1-alpha.seq[k]/2)
      ci.l.k = gx.xseq.b[,x0] - quant.k*se.xseq.b[,x0]
      ci.r.k = gx.xseq.b[,x0] + quant.k*se.xseq.b[,x0]
      pi.hat[k] = mean(1*(gx.xseq[x0]<=ci.r.k & gx.xseq[x0]>=ci.l.k))
    }
    ind.min = which.min(abs(pi.hat-cov))
    alpha.hat[x0] = alpha.seq[ind.min]
  }
  
  qz.hh = qnorm(1-quantile(alpha.hat,1-vareps))
  out = list(gx.hh=gx.hh, se.hh=se.hh, qz.hh=qz.hh)
}


W.fun = function(u,kernel){
  if (kernel=="epanechnikov" | kernel=="epa") w = 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel=="uniform"      | kernel=="uni") w =          0.5*(abs(u)<=1)
  if (kernel=="triangular"   | kernel=="tri") w =   (1-abs(u))*(abs(u)<=1)
  return(w)  
}

k.fun = function(u){
  if (kernel=="epanechnikov" | kernel=="epa") w = 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel=="uniform"      | kernel=="uni") w =          0.5*(abs(u)<=1)
  if (kernel=="triangular"   | kernel=="tri") w =   (1-abs(u))*(abs(u)<=1)
  return(w)  
}

qrXXinv = function(x, ...) {
  #tcrossprod(solve(qr.R(qr(x, tol = 1e-10)), tol = 1e-10))
  #tcrossprod(solve(qr.R(qr(x))))
  chol2inv(chol(crossprod(x)))
}
