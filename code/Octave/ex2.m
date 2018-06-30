clear all;
source tpunc.m;
msg = @(s) fprintf("%s\n",s);
rep = @(s,d) fprintf("%s: %d\n",s,d);
more off;

M = 50;
t = 0.1:0.1:4.0;

ESV(1:length(t))=0;
EKV(1:length(t))=0;

cc=0;

for R=t
  H=@(x,k) int8(x^2+k^2>R);
  cc=cc+1;
  rep("Step",cc);

  # Find the Schroedinger ground state for H
  [phi0, ESV(cc), info, iter, nf, lambda,x] = findGroundState(H,M);

  # Ensure boundary condition
  phi0(1)=0;phi0(end)=0;

  # Find Kantorovich min for ground state phi0 using 'glpk'
  [gamma1,mu1,nu1,x,k,EKV(cc),ES1,HM,phi1,FMIN,ERRNUM] = SolveWithGLPK(H,phi0',M);

  # Find Kantorovich min for ground state phi0 using 'linprog'
  [gamma2,mu2,nu2,x,k,EKV2(cc),ES2,HM,phi2,FMIN] = SolveWithLinprog(H,phi0',M);

  # Check results
  #msg("Check results ...");
  #checkGLPKSolution(gamma1,mu1,nu1,EKV(cc),ERRNUM,FMIN);
  #checkGLPKSolution(gamma2,mu2,nu2,EK2,0,EK2);
  clear phi0 info iter nf lambda x gamma1 mu1 nu1 k ES1 HM phi FMIN ERRNUM
endfor

save data4 t ESV EKV EKV2

# plot(t,ESV,t,EKV)






