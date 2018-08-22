source tpunc.m;

function [phi0, obj, info, iter, nf, lambda,x] = GS(H,M,f0)
  [N,L,x,k,dx,dk] = setup_grid(M);
  HM = FunctionMatrix(H,x,k); 
  U = UCDFT(M);
  E = @(phi) ESchroed(HM,U,M,phi,x,k); 
  X0 = makePhi(f0,x);
  [phi0, obj, info, iter, nf, lambda] = sqp (X0, E, @g, []);
endfunction


#f0 = @(t) int8(abs(t)<1);
f0 = @(t) 1/(1+t^2); 
