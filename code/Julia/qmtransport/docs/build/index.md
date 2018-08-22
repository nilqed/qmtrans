
<a id='qmt.jl-Documentation-1'></a>

# qmt.jl Documentation



<a id='qmt.Grid' href='#qmt.Grid'>#</a>
**`qmt.Grid`** &mdash; *Type*.



```
	Grid(N,L,dx,dk,x,k)
```

The Grid `type` is a `struct` that stores the values of a two-dimensional grid, where `N::Int` denotes the number of points, `L::Real` the length, `dx::Real` and `dk::Real` the spacing and finally `x::Real,k::Real` the vectors of points. Although not checked, it is expected that `L=N*dx=N*dk` and `length(x)=N`.

<a id='qmt.setupGrid' href='#qmt.setupGrid'>#</a>
**`qmt.setupGrid`** &mdash; *Function*.



```
	setupGrid(M::Int)
```

Given a positive integer `M`, return a `Grid` with the following values:

```
	N = 2*M+1
  L = sqrt(2*N*pi)
	dx = dk=L/N
	x = dx * [j for j in -M:M]
	k = dk * [j for j in -M:M]
```

