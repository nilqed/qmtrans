module qmt
export Grid, setupGrid
#=
 Description:
 eee
 eee
=#

println("\nJulia -- qmtransport, grid\n\n")
using Plots

"""
		Grid(N,L,dx,dk,x,k)

The Grid `type` is a `struct` that stores the values of a two-dimensional grid,
where `N::Int` denotes the number of points, `L::Real` the length, `dx::Real`
and `dk::Real` the spacing and finally `x::Real,k::Real` the vectors of points.
Although not checked, it is expected that `L=N*dx=N*dk` and `length(x)=N`.
"""
struct Grid
    N::Int
    L::Real
    dx::Real
    dk::Real
    x::Vector{Real}
    k::Vector{Real}
end

"""
		setupGrid(M::Int)

Given a positive integer `M`, return a `Grid` with the following values:
```
	N = 2*M+1
  L = sqrt(2*N*pi)
	dx = dk=L/N
	x = dx * [j for j in -M:M]
	k = dk * [j for j in -M:M]
```
"""
function setupGrid(M::Int)::Grid
  N=2*M+1
  L=sqrt(2*N*pi)
  dx=L/N
  dk=dx
  x=dx * [j for j in -M:M]
  k=dk * [j for j in -M:M]
  return Grid(N,L,dx,dk,x,k)
end

"""
		deMoivreNumber(g:Grid)::Complex

Compute the standard *de Moivre number* of the grid `g`, that is
```
	exp(2*pi*im/g.N) = cos(2*pi/g.N) + im * sin(2*pi/g.N)
```
See [WP:Root of unity](https://en.wikipedia.org/wiki/Root_of_unity).
"""
function deMoivreNumber(g::Grid)::Complex
  return cos(2*pi/g.N) + im * sin(2*pi/g.N)
end

"""
		hamiltonMatrix(H::Function,g::Grid)::Matrix{Real}

Given a classical Hamilton function `H(x,k)` on the grid `g` (phase space),
compute the corresponding matrix comprising the values of `H` at the grid
points `(g.x,g.k)`.
"""
function hamiltonMatrix(H::Function,g::Grid)::Matrix{Real}
  return [H(x,k) for x in g.x, k in g.k]
end

"""
		func2vector(f::Function,g::Grid)::Vector{Complex}

Given a (complex valued) function on the real line, compute the corresponding
vector comprising the values at the grid points `g.x`. Here, `f` usually is
a Schroedinger wave function.
"""
function func2vector(f::Function,g::Grid)::Vector{Complex}
  return [f(x) for x in g.x]
end

"""
		signum(n::Int)::Int

Compute the `signum` for the integer `n`, that is
```
	signum(n) = (-1)^n.
```
"""
function signum(i::Int)::Int
  if isodd(i)
    return(-1)
  else
    return(1)
  end
end

"""
		ucdftMatrix(g::Grid)::Matrix{Complex}

Compute the unitary centered discrete Fourier transform matrix for the
grid `g`. Note: the inverse transform is given by the conjugate matrix
(use the conjugate transposition operator `'` or simply  `conj`
by symmetry).
"""
function ucdftMatrix(g::Grid)::Matrix{Complex}
	if ~isodd(g.N) error("Grid.N must be odd!") end
	M =(g.N-1)/2
	a = deMoivreNumber(g)
  F = [a^((m-M-1)*(n-M-1)) for n in 1:g.N, m in 1:g.N]
  return 1/sqrt(g.N) * F
end


"""
		stdHamiltonian(g::Grid,Vf::Function,h::Real,m::Real)::Matrix{Real}

Compute the standard Hamiltonian with potential function `Vf`, mass `m` and
Planck's constant `h` (value depends on the unit system chosen). The
eigenvalues and eigenfunctions of the resulting matrix represent the
quantum mechanical states of the system `p^2/(2*m)+Vf(x)`. For details
see
>[1] Gabriel G. Balint-Kurti, Richard N. Dixon, and C. Clay Marston. Grid
>methods for solving the schr¨ odinger equation and time dependent quan-
>tum dynamics of molecular photofragmentation and reactive scattering
>processes. International Reviews in Physical Chemistry, 11(2):317–344, 1992.
and
>[2] Kurt Pagani, Transference Plans and Uncertainty,
>arXiv,[math.AP]
"""
function stdHamiltonian(g::Grid,Vf::Function,h::Real,m::Real)::Matrix{Real}
  N = g.N; L = g.L
  f = (i,j) -> (signum(i-j)/m * cos(pi*(i-j)/N) * (h/(2*L*sin(pi*(i-j)/N)))^2)
  k = h^2/(4*m*L^2)*(N^2-1)/6
  V=func2vector(Vf,g)
  T=eye(N)
  [T[i,j]=f(i,j) for i=1:N, j=1:N if i!=j]
  [T[i,i]=k+V[i] for i=1:N]
  return T
end

"""
	checkEigResult(M,eigval,eigvec)

Check the quality of the solution to the eigenvalue problem
```
  M * eigvec[1:n,i] = eigval[i] * eigvec[1:n,i]
```
where `eigval=eigvals(M)` and `eigvec=eigvecs(M)`.
The return value is an array of the norms of `(M-e)v`, so that
```
  err = sum(checkEigResult(M,ev,EV))
```
should be near to zero.
"""
function checkEigResult(M,eigval,eigvec)
  n=size(M)[1]
  E=eye(n)
  return [norm((M-eigval[i]*E)*eigvec[1:n,i]) for i=1:n]
end


"""
		QMDATA(g,h,m,H,ev,EV,err)

The type `QMDATA` stores the results of an eigenvalue/eigenvector calculation
on the grid `g`.
"""
struct QMDATA
		g::Grid
		h::Real
		m::Real
		H::Matrix{Real}
		ev::Vector{Real}
		EV::Matrix{Real}
		err::Vector{Real}
end

"""
		solveSchroedinger(M::Int,V::Function,h::Real,m::Real)::QMDATA

Solve the (time independent) Schroedinger equation on a N=2*M+1 grid, where
`V` is the potential function, `m` the mass and `h` the Planck constant in
suitable units. The return value is a `QMDATA` struct.
"""
function solveSchroedinger(M::Int,V::Function,h::Real,m::Real)::QMDATA
	g = setupGrid(M)
	H = stdHamiltonian(g,V,h,m)
	ev = eigvals(H)
	EV = eigvecs(H)
	err = checkEigResult(H,ev,EV)
	return QMDATA(g,h,m,H,ev,EV,err)
end

"""
		exHarmonicOscillator(M::Int)::QMDATA

This is an example for `solveSchroedinger` with `V(x)=x^2` and
`h=2*sqrt(2)*pi` and `m=1`, such that `hbar^2/(2*m)=1`. Therefore, the
eigenvalues and eigenvectors are those of the harmonic oscillator and
exactly known:
```
  E_n = 2*n+1, n=0,1,2,...
	psi_n(x) = hermite(n,x)
```
"""
function exHarmonicOscillator(M::Int)::QMDATA
	h = 2*sqrt(2)*pi
	m = 1
  return solveSchroedinger(M,x->x^2,h,m)
end

"""
		plotEigenvector(q,i,pm=gr)

Plot eigenvector `i` from q::QDATA using the method
```
   pr = gr | plotly | plotlyjs
```
Note that in order to plot in the Juno pane, `plotlyjs` is required, whereas
`plotly` will plot to the browser (e.g. Firefox). Probably you need to
install the package: `Pkg.add("PlotlyJS")`
"""
function plotEigenvector(q::QMDATA,i::Int,pm::Function=gr)
	pm()
	p = plot(q.g.x,q.EV[1:q.g.N,i]);
	return p
end


function hermitePolynmial(n::Int,x::Real)::Real
	fac = k->factorial(k)
	s=[signum(m)/fac(m)/fac(n-2*m)*(2*x)^(n-2*m) for m in 0:Int(floor(n/2))]
	return factorial(n)*sum(s)
end

function hermite(n::Int,x::Real)::Real
	c = 1/sqrt(2^n*factorial(n)*sqrt(pi))
	return c*hermitePolynmial(n,x)*exp(-x^2/2)
end

"""
		schroedingerEnergy(H::Function,phi::Function,g::Grid)::Real

Calculate the Schroedinger energy for the Hamilton function `H` and the
wave function `phi` on the grid `g`. That is
```
	E(phi) = SUM {H(x,k) |phi(x)|^2 |UCDFT(phi)(k)|^2 : x,k in g}
```
# Examples
```julia-repl
julia> g = setupGrid(50);
julia> schroedingerEnergy((x,k)->x^2+k^2,x->hermite(0,x),g)
1.0000000000000009
julia> schroedingerEnergy((x,k)->x^2*k^2,x->hermite(0,x),g)
0.2500000000000002
```
"""
function schroedingerEnergy(H::Function,phi::Function,g::Grid)::Real
	cdot = x -> x*conj(x)
  HM = hamiltonMatrix(H,g)
	phiv = func2vector(phi,g)
	phiv = phiv/norm(phiv) # normalized => sum(phiv)=1
	v = map(cdot,phiv)
	ftv = map(cdot,ucdftMatrix(g)*phiv) # sum(ftv)=1 ?
	return Real(dot(ftv,HM*v))
end


# R = exHarmonicOscillator(50)
# sum(R.err)
# p1 = plotEigenvector(R,1); p1
# plotEigenvector(R,1,gr)
# plotEigenvector(R,1,plotly)
# plotEigenvector(R,1,plotlyjs)


### Example (tested with M=100 => N=201)

#=
g=setupGrid(5)
HM=hamiltonMatrix((x,k)->x^2+k^2,g)
size(HM)

h=2*sqrt(2)*pi
H0=stdHamiltonian(g,x->x^2,h,1)

ev=eigvals(H0)
EV=eigvecs(H0)
err = sum(checkEigResult(H0,ev,EV))
=#



#using Plots
#gr()
#plot(g.x,ev)

#=
using Plots
plotly()
plot(g.x,EV[1:g.N,8])
=#

#l=(x,k)->x^2+k^2

#makeGrid(N::Int,dx::Real)::Grid = Grid(N,dx)

#makeGrid(12,0.002)
# --> Grid(12, 0.002)
#ans.N
# --> 12

#h=1;
#m=1;
#L=1;
#N=10;

#f(i,j) = ((-1)^(i-j)/m * (h/(2*L*sin(pi*(i-j)/N)))^2);
#g = h^2/(4*m*L^2)*((N-1)*(N-2)/6+(N/2));

#T=[(0.0;0.0) for i=1:N, j=1:N];
#[T[i,j]=f(i,j) for i=1:N, j=1:N if i!=j];
#[T[i,i]=g for i=1:N];

# psi0(x) = pi^(-1/4) * exp(-x^2/2) -- Hermite0
# P = func2vector(psi0,g)
# plot(g.x,[-EV[1:g.N,1],EV[1:g.N,2],EV[1:g.N,3],real(P)])
# plot(g.x,[-EV[1:g.N,1],EV[1:g.N,2],EV[1:g.N,3],real(P)/norm(P)])
# plot(g.x,[-EV[1:g.N,1],real(P)/norm(P)])


#quit()

end