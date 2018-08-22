library(transport)

setClass("grid", slots = 
  list(N  = "numeric",
       L  = "numeric",
      dx  = "numeric",
      dk  = "numeric",
       x  = "numeric",
       k  = "numeric", 
       deMoivreNumber = "complex",
       ucdftMatrix = "matrix"))

#grid = new("Grid", N = 10)

grid = function(M) {
  N=2*M+1
  L = sqrt(2*N*pi)
  dx = L/N
  dk = dx
  x = dx * (-M:M)
  k = dk * (-M:M)
  dmn = exp(2*pi*1i/N)
  sqn = 1/sqrt(N)
  dft = matrix(0,N,N)
  for (i in 1:N) {for (j in i:N)
  {dft[i,j]=dft[j,i]=sqn*dmn^((i-M-1)*(j-M-1))}}
  new("grid", 
      N = as.integer(N),L=as.numeric(L),
      dx = as.numeric(dx), dk=as.numeric(dk),
      x = as.numeric(x), k = as.numeric(k),
      deMoivreNumber = as.complex(dmn),
      ucdftMatrix = as.matrix(dft))
}




#g=grid(100)
#g@deMoivreNumber^201  # we have @ for S4 and $ for S3
# list comprehensions in R? Yes, why not
# l=1:10, l[l>4]
# l[l==1 | l==4]
# U=.Last.value@ucdftMatrix or
# U=g@ucdftMatrix
# ev=eigen(U)$values
# EV=eigen(U)$vectors
# prod(ev)
# abs(prod(ev)) # should be 1
# Conj(U) # complex conjugate


is.even = function(x) x %% 2 == 0 
is.odd = function(x) x %% 2 != 0 
signum = function(x) if (is.odd(x)) {-1} else {1}




setClass("schroed", slots =
  list(V = "numeric",
       h = "numeric",
       m = "numeric",
       H = "matrix",
       ev = "numeric",
       EV = "matrix"),
  contains = "grid")

#  typeof(s) -> S4
#  class(s) -> "schroed" attr(,"package")
#  s=schroed(5)
#  > s@N
#  [1] 11
# t=schroed(2,sin)

schroed = function(M,V,h=2*sqrt(2)*pi,m=1) {
  g = grid(M)
  N = g@N
  if (is.function(V)) {vecV=V(g@x)} else {vecV=V} #more checks
  H = matrix(0,N,N)
  k = h^2/(4*m*g@L^2)*(N^2-1)/6
  for (i in 1:N) 
    {for (j in i:N)
      { if (j != i) 
          {a = (i-j)*pi/N
           H[i,j]=H[j,i]=signum(i-j)/m*cos(a)*(h/(2*g@L*sin(a)))^2}
        else
        {H[i,i]=k+vecV[i]}
    }
  }
  eig = eigen(H)
  ev = eig$values  # S3 objects
  EV = eig$vectors
  new("schroed",g,V=vecV,h=h,m=m,H=H,ev=ev,EV=EV)
}

# hosc=schroed(20,function(x){x^2})
# H=hosc@H
# ev=eigen(H)@values
# EV=eigen(H)$vectors
# plot(hosc@x,-EV[,41])
# plot(hosc@x,-EV[,41],type="l")
# plot(hosc@x,-EV[,41],type="h")
# ...
# hosc@ev
# > length(hosc@H) for M=200
# [1] 160801> length(hosc@H)

setGeneric("displayEV",
  function(obj,n)
  standardGeneric("displayEV"))

setMethod("displayEV", signature(obj = "schroed", n = "numeric"),
  function(obj,n)
    plot(obj@x,obj@EV[,n],type="l")
  )

#> displayEV(hosc,1)
#> displayEV(hosc,hosc@N)
# p=exp(-hosc@x^2/2)/(pi^(1/4))
# np=sqrt(sum(p*p))
# phi0=p/np
# plot(hosc@x,phi0-hosc@EV[,hosc@N],type="l")
#> d=phi0-hosc@EV[,hosc@N]
#> dd=sqrt(sum(d*d))
#> dd
#[1] 8.935168e-14
# eigenvalues: plot(1:hosc@N,rev(hosc@ev))
# plot(1:hosc@N,rev(hosc@ev),type="h")
# plot(1:30,rev(hosc@ev)[1:30],type="h")

# V=function(x){as.integer(abs(x) < 4)*(-10)}
# pot=schroed(100,V)
# displayEV(pot,201)
# displayEV(pot,197)
# pot@ev -- bound states?


setGeneric("checkEV",
           function(obj)
             standardGeneric("checkEV"))

setMethod("checkEV", signature(obj = "schroed"),
          function(obj){
            s=0
            for (n in 1:obj@N)
            {s=s+norm(obj@H%*%obj@EV[,n]-obj@ev[n]*obj@EV[,n])}
            return(s)}
)

# checkEV(hosc) 

hamiltonMatrix0 = function(H,g) {
  n = g@N
  M = matrix(0,n,n)
  for (i in 1:n) {
    for (j in 1:n) {
      M[i,j]=H(g@x[i],g@k[j])
    }
  }
  return(M)}

#> HM=hamiltonMatrix0(H,g) -- better see below

hamiltonMatrix = function(H,g) {outer(g@x,g@k,FUN=H)}
  
func2vector = function(f,g) {f(g@x)}

#> f0=function(x){sin(x)}
#> v=func2vector(f0,g)
#> plot(g@x,v,type="l")

# HM=hamiltonMatrix(H,g)
# persp(g@x,g@k,HM)


setClass("kantorovich", slots =
           list(H = "function",
                phi = "function",
                mu = "numeric",
                nu = "numeric",
                HM = "matrix",
                EK = "numeric"),
         contains = "grid")

kantorovich = function(M,H,phi) {
  g = grid(M)
  N = g@N
  if (is.function(phi)) {v=phi(g@x)} else {v=as.vector(phi)}
  if (is.function(H)) {HM=outer(g@x,g@k,FUN=H)} else {HM=as.matrix(H)}
  v = v/sqrt(sum(v*Conj(v)))
  w = g@ucdftMatrix %*% v
  mu = v*Conj(v)
  nu = w*Conj(w)
  new("kantorovich",g,H=H,phi=phi,mu=as.numeric(mu),nu=as.numeric(nu),
      HM=as.matrix(HM))
  }

# NOTE: grid g must be added in 'new' if inheritance is wished
# H=function(x,k){x^2+k^2}
# K=kantorovich(20,H,function(x){exp(-x^2)})
#> plot(K@x,K@mu,type="l")
#> plot(K@k,K@nu,type="l")
# EK=K@nu %*% (K@HM %*% K@mu)  #TODO
# 
#> pm=K@mu %*% t(K@nu)  ## product measure
#> persp(K@x,K@k,pm)    ## plot it 3d
#> contour(K@x,K@k,pm)  ## contour
#> contour(K@x,K@k,pm,nlevels=8,xlim=range(-1,1),ylim=range(-3,3))
#> contour(K@x,K@k,pm,nlevels=4,xlim=range(-3,3),ylim=range(-3,3))

K=kantorovich(60,H,function(x){exp(-x^2/2)})
T=transport(K@mu,K@nu,K@HM,method="revsimplex")

gammaM = function(N,T) {
  gm=matrix(0,N,N)
  mm=length(T$from)
  for (i in 1:mm) {
    gm[T$from[i],T$to[i]]=T$mass[i]
  }
  return(gm)
}
   
gamma=gammaM(K@N,T)      
# length(gamma[gamma>0]) # n+m-1 
# contour(K@x,K@k,gamma)
#> contour(K@x,K@k,gamma,nlevels=8,xlim=range(-1,1),ylim=range(-3,3))
#> contour(K@x,K@k,gamma,nlevels=8,xlim=range(-2,2),ylim=range(-3,3))
#> sum(apply(gamma, 1, sum))
#[1] 1
#> sum(apply(gamma, 2, sum))
#[1] 1
# sum(abs(apply(gamma, 1, sum)-K@mu))

#> sum(diag(t(K@HM) %*% gamma))
#[1] 1

transp <- matrix(0,K@N,K@N)
transp[cbind(T$from,T$to)] <- T$mass
rownames(transp) <- paste(ifelse(nchar(K@mu)==2," ",""),K@mu,sep="")
colnames(transp) <- paste(ifelse(nchar(K@nu)==2," ",""),K@nu,sep="")
#print(transp)

#contour(K@x,K@k,transp,nlevels=12,xlim=range(-2,2),ylim=range(-3,3))
#sum(diag(t(K@HM) %*% t(transp))) ---> 1
