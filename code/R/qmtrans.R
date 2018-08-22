#####################################################################
# S4 Class 'psgrid' (typeof -> S4)
# (2*M+1)*(2*M+1) point grid, symmetric to (0,0).
# Construction: g=psgrid(M)
#####################################################################

setClass("psgrid", slots = 
           list(N  = "numeric",
                L  = "numeric",
                dx  = "numeric",
                dk  = "numeric",
                x  = "numeric",
                k  = "numeric",
                deMoivreNumber = "complex",
                ucdftMatrix = "matrix",
                dftv   = "function",
                dftf   = "function",
                idftv  = "function",
                idftf  = "function",
                eval   = "function",
                evplot = "function",
                plotv  = "function",
                show   = "function"))

###
### Constructor psgrid(PositiveInteger)
###
psgrid = function(M) {
  N=2*M+1
  L = sqrt(2*N*pi)
  dx = L/N
  dk = dx
  x = dx * (-M:M)
  k = dk * (-M:M)
  dmn = exp(2*pi*1i/N)
  sqn = 1/sqrt(N)
  #
  dft = matrix(0,N,N)
  for (i in 1:N) {for (j in i:N)
  {dft[i,j]=dft[j,i]=sqn*dmn^((i-M-1)*(j-M-1))}}
  xeval = function(f){f(x)}
  keval = function(f){f(k)}
  #
  evplot = function(f,t='l') {v=xeval(f); plot(x,v,type=t);v}
  plotv = function(v,t='l'){plot(x,v,type=t)}
  dftv = function(v) {dft %*% v}
  dftf = function(f) {dft %*% xeval(f)}
  idftv = function(v) {solve(dft) %*% v}
  idftf = function(f) {solve(dft) %*% keval(f)}
  #
  show = function(colour=1){
    plot(x,k,type='n',xlab="position-x",ylab="mommentum-k")
    graphics::grid(N,N,col=colour)}
  #
  new("psgrid", 
      N = as.integer(N),L=as.numeric(L),
      dx = as.numeric(dx), dk=as.numeric(dk),
      x = as.numeric(x), k = as.numeric(k),
      deMoivreNumber = as.complex(dmn),
      ucdftMatrix = as.matrix(dft),
      dftv = as.function(dftv),
      dftf = as.function(dftf),
      idftv = as.function(idftv),
      idftf = as.function(idftf),      
      eval = as.function(xeval),
      evplot = as.function(evplot),
      plotv = as.function(plotv),
      show = as.function(show))
}

###
### Methods
###

setGeneric("showGrid",
           function(obj,col)
             standardGeneric("showGrid"))

setMethod("showGrid", signature(obj = "psgrid", col = "numeric"),
          function(obj,col){
            plot(obj@x,obj@k,type='n',xlab="position-x",ylab="mommentum-k")
            graphics::grid(obj@N,obj@N,col=col)})


setGeneric("func2vector",
           function(obj,f)
             standardGeneric("func2vector"))

setMethod("func2vector", signature(obj="psgrid", f="function"),
          function(obj,f){f(obj@x)})


setGeneric("selfTest",
           function(obj)
             standardGeneric("selfTest"))

setMethod("selfTest", signature(obj="psgrid"),
          function(obj){
            obj@show()
            A=obj@dftf(sin)
            B=obj@idftv(A)
            C=obj@eval(sin)
            sum(abs(C-Re(B)))<0.0000001
          })


#####

is.even = function(x) x %% 2 == 0 
is.odd = function(x) x %% 2 != 0 
signum = function(x) if (is.odd(x)) {-1} else {1}


#####################################################################
# S4 class "schroed1d" 
# Solve the 1d Schroedinger equation -Lap(f)+Vf=Ef
# Constructor: s=schroed1d(M:PositiveInteger,V:Function)
#####################################################################

setClass("schroed1d", slots =
           list(V = "numeric",
                h = "numeric",
                m = "numeric",
                H = "matrix",
                ev = "numeric",
                EV = "matrix",
                plotEigVec = "function"),
         contains = "psgrid")

###
### Constructor schroed1d(M:PositiveInteger,V:Function,h,m) 
###

schroed1d = function(M,V,h=2*sqrt(2)*pi,m=1) {
  g = psgrid(M)
  N = g@N
  if (is.function(V)) {vecV=V(g@x)} else {vecV=V} #more checks
  H = matrix(0,N,N)
  k = h^2/(4*m*g@L^2)*(N^2-1)/6
  #
  for (i in 1:N) 
  {for (j in i:N)
  { if (j != i) 
  {a = (i-j)*pi/N
  H[i,j]=H[j,i]=signum(i-j)/m*cos(a)*(h/(2*g@L*sin(a)))^2}
    else
    {H[i,i]=k+vecV[i]}}}
  #
  eig = eigen(H)
  ev = eig$values  # S3 objects
  EV = eig$vectors
  pEV = function(n,t='l',mt=paste("N=",toString(N),"dx=",
          toString(round(g@dx,3)))) 
  {plot(g@x,EV[,n],type=t,
        xlab="x",ylab=paste("Eigenfunction#",toString(n)))
    title(mt,paste("eigval=",toString(round(ev[n],4))))}
  #
  new("schroed1d",g,V=vecV,h=h,m=m,H=H,ev=ev,EV=EV,plotEigVec=pEV)
}

###
### Methods
###

setGeneric("displayEV",
           function(obj,n)
             standardGeneric("displayEV"))

setMethod("displayEV", signature(obj = "schroed1d", n = "numeric"),
          function(obj,n)
            plot(obj@x,obj@EV[,n],type="l"))

setGeneric("checkEV",
           function(obj)
             standardGeneric("checkEV"))

setMethod("checkEV", signature(obj = "schroed1d"),
          function(obj){
            s=0
            for (n in 1:obj@N)
            {s=s+norm(obj@H%*%obj@EV[,n]-obj@ev[n]*obj@EV[,n])}
            return(s)})



#####################################################################
# S4 class "kantorovich" 
# Calculate the Kantorovich energy inf{Trace(H*g): g in Gamma(phi)}
# Constructor: k=kantorovich(M:PositiveInteger,H:func,phi:func)
#####################################################################

library(transport)

setClass("kantorovich", slots =
           list(H = "function",
                phi = "function",
                mu = "numeric",
                nu = "numeric",
                HM = "matrix",
                ES = "numeric",
                EK = "numeric",
                prodMeasure = "matrix",
                transfPlan = "list",
                gamma = "matrix",
                display = "function"),
         contains = "psgrid")


kantorovich = function(M,H,phi) {
  g = psgrid(M)
  N = g@N
  if (is.function(phi)) {v=phi(g@x)} else {v=as.vector(phi)}
  if (is.function(H)) {HM=outer(g@x,g@k,FUN=H)} else {HM=as.matrix(H)}
  v = v/sqrt(sum(v*Conj(v)))
  w = g@dftv(v)
  mu = as.double(v*Conj(v))
  nu = as.double(w*Conj(w))
  ES = as.double(sum(nu * (HM %*% mu)))
  pm = as.matrix(mu %*% t(nu))
  T=transport(mu,nu,HM,method="revsimplex")
  gammaM = function(N,T) {
    gm=matrix(0,N,N)
    mm=length(T$from)
    for (i in 1:mm) {
      gm[T$from[i],T$to[i]]=T$mass[i]
    }
    return(gm)}
  gamma=as.matrix(gammaM(N,T))
  EK = sum(diag(HM%*%t(gamma))) #tr(H*g')
  display = function(scale=10) {
    r=g@L/scale
    rx=range(-r,r)
    ry=range(-r,r)
    par(mfrow=c(2,2))
    ylx=expression(paste(mu,"(x)"))
    ylk=expression(paste(nu,"(k)"))
    t1=expression(paste(mu," x ",nu))
    t2=expression(paste(gamma))              
    plot(g@x,mu,type="h",col="green",xlab="x",ylab=ylx)
    plot(g@k,nu,type="h",col="red",xlab="k",ylab=ylk)
    contour(g@x,g@k,pm,xlim=rx,ylim=ry,col="blue",main=t1)
    contour(g@x,g@k,gamma,xlim=rx,ylim=ry,col="brown",main=t2)}
  # create
  new("kantorovich",g,H=H,phi=phi,mu=mu,nu=nu,
      HM=as.matrix(HM),ES=ES,prodMeasure=pm, transfPlan=T,
      gamma=gamma,EK=as.numeric(EK),display=as.function(display))
}

###
