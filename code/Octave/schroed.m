clear all;

#H = @(x,k) ((x^2+k^2>10));  # 5,10,20 -> ???

#H = @(x,k) (x^2+k^2);  # M=50 almost perfect match
#H = @(x,k) (x^2*k^2);  # M=50 almost perfect match

delta=@(x) 1000*(x==0);

function pot = V(x)
  if (abs(x)<10)
    pot=60;  
  else
    pot=0;
  endif
endfunction

H = @(x,k) (k^2+V(x));

M = 50

N = 2*M+1
L = sqrt(2*pi*N)
dx = L/N
dk = 2*pi/L

#x = linspace(-M,M,N) 
x = (-M:1:M)*dx;
k = (-M:1:M)*dk;


# MkMatrix
function retval = MkMatrix(f,x,y)
  # MkMatrix(f,x,y)= matrix f(x_i,y_j)
  N = size(x)(2);
  F = eye(N);
  for m=1:N
    for n=1:N
      F(m,n) = f(x(m),y(n));
    endfor
  endfor
  retval = F;
endfunction

# Unitary Centered Discrete Fourier Transform (UCDFT)
function retval = UCDFT(M)
  # M = (N-1)/2 : integer
  # Computes the NxN matrix 
  #   UCDFT(M)= 1/sqrt(pi)*[exp(2*i*pi/N)^((m-M-1)*(n-M-1))]
  # 
  N = 2*M+1;
  U = eye(N);
  w = exp(-2*i*pi/N);
  for m=1:N
    for n=1:N
      U(m,n)=w^((m-M-1)*(n-M-1));
    endfor
  endfor
  retval = U/sqrt(N);
endfunction


# Compute the UCDFT matrix
global U = UCDFT(M);
global HM = MkMatrix(H,x,k);


# Energy function
function obj = E(phi)
  global U HM;
  X = phi'/norm(phi);
  mu = abs(X).^2;
  nu = abs(U*X').^2;
  g  = mu.*nu;
  obj = trace(HM*g');
endfunction

# Conditions: X(1)=0, X(N)=0, |X|=1.
function constraint = g(phi)
  a = size(phi)(2);
  b = size(phi)(1);
  constraint = [sumsq(phi)-1; phi(a); phi(b)];
endfunction


# Test 1
X(1:N)=1;
X(1)=X(N)=0;
E(X')

[phi0, obj, info, iter, nf, lambda] = sqp (X, @E, @g, []);
E(phi0)
h0 = exp(-x.^2/2)/norm(exp(-x.^2/2));
plot(x,phi0/norm(phi0),x,h0);

  