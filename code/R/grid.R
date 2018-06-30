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
