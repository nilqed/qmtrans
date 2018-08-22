###############################################################################
# 
# 
#  
#

clear all;

M = 100;

function [N,L,x,k,dx,dk] = setupOddGrid(M)
  # Set up a grid [-M,..,-1,0,1,...,M] with N=2*M+1 (odd number of grid points)
  # The length L is determined by the requirement dx=dk, i.e L^2=2*pi*N.
  N = 2*M+1;
  L = sqrt(2*pi*N);
  dx = L/N;
  dk = 2*pi/L;
  x = (-M:1:M)*dx; #x = linspace(-M,M,N)*dx 
  k = (-M:1:M)*dk; #k = linspace(-M,M,N)*dk 
endfunction

function [N,L,x,k,dx,dk] = setupEvenGrid(M)
  # Set up a grid [-M+1,..,-1,0,1,...,M] with N=2*M (even number of grid points)
  # The length L is determined by the requirement dx=dk, i.e L^2=2*pi*N.
  N = 2*M;
  L = sqrt(2*pi*N);
  dx = L/N;
  dk = 2*pi/L;
  x = (-(M-1):1:M)*dx; #x = linspace(-M,M,N)*dx 
  k = (-(M-1):1:M)*dk; #k = linspace(-M,M,N)*dk 
endfunction

function [HM,phi,E] = fghodd(M,V,m,h)
  # Solve the Schroedinger equation by the FGH method
  # M is a positive integer (determines N=2*M+1=#grid points)
  # V is the potential function 
  # m is the mass of the particle
  # h is the Planck constant (in suitable units)
  # Returns: 
  #  HM is the Hamiltonian matrix, 
  #  phi are the eigenvectors
  #  E are the eigenvalues (sorted)  
  [N,L,x,k,dx,dk] = setupOddGrid(M);
  VV = arrayfun(V,x);
  fac = h^2/(8*pi*m*N);
  HM = eye(N);
  for m=1:N
    for n=1:N
      if (m==n)
        HM(m,n) = fac*(N^2-1)/6 + VV(m);
      else
        HM(m,n) = fac*(-1)^(m-n)*cos(pi*(m-n)/N)/sin(pi*(m-n)/N)^2;
      endif
    endfor
  endfor
  [vv,eed] = eig(HM);
  [E,perm] = sort(diag(eed));
  phi = vv(:,perm);
endfunction


### Examples
function [HM,phi,E] = harmonic_oscilator(M)
  [N,L,x,k,dx,dk] = setupOddGrid(M);
   V = @(t) t^2;
  [HM,phi,E] = fghodd(M,V,1/2,2*pi);  # hbar^2=2*m => m=1/2, h=2*pi
   g0 = @(t) exp(-t^2/2);
   h0 = arrayfun(g0,x);
   h0 = h0/norm(h0);
   plot(x,-phi(:,1),x,h0,x,phi(:,2));
endfunction

function [HM,phi,E] = potential_well(M)
  [N,L,x,k,dx,dk] = setupOddGrid(M);
   V = @(t) -5*int8(abs(t)<4);
  [HM,phi,E] = fghodd(M,V,1/2,2*pi);  # hbar^2=2*m => m=1/2, h=2*pi
   g0 = @(t) exp(-t^2/2);
   h0 = arrayfun(g0,x);
   h0 = h0/norm(h0);
   plot(x,phi(:,1),x,phi(:,2),x,phi(:,3));
endfunction

function [HM,phi,E] = box(M)
  [N,L,x,k,dx,dk] = setupOddGrid(M);
   V = @(t) 1000*int8(abs(t)>4);
  [HM,phi,E] = fghodd(M,V,1/2,2*pi);  # hbar^2=2*m => m=1/2, h=2*pi
   g0 = @(t) exp(-t^2/2);
   h0 = arrayfun(g0,x);
   h0 = h0/norm(h0);
   plot(x,phi(:,1),x,phi(:,2),x,phi(:,3));
endfunction

function [HM,phi,E] = step_potential(M)
  [N,L,x,k,dx,dk] = setupOddGrid(M);
   V = @(t) 10*int8(t>=0);
  [HM,phi,E] = fghodd(M,V,1/2,2*pi);  # hbar^2=2*m => m=1/2, h=2*pi
   g0 = @(t) exp(-t^2/2);
   h0 = arrayfun(g0,x);
   h0 = h0/norm(h0);
   plot(x,phi(:,1),x,phi(:,2),x,phi(:,3));
endfunction


