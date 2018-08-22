clear all;
source tpunc.m;
msg = @(s) fprintf("%s\n",s);
rep = @(s,d) fprintf("%s: %d\n",s,d);
more off;

M =100;
H=@(x,k) int8(x^2+k^2>1);

# Find the Schroedinger ground state for H
msg("Find ground state ...");
[phi0, ES0, info, iter, nf, lambda,x] = findGroundState(H,M);

# Ensure boundary condition
phi0(1)=0;phi0(end)=0;

# Find Kantorovich min for ground state phi0 using 'glpk'
msg("Find gamma1 ...");
[gamma1,mu1,nu1,x,k,EK1,ES1,HM,phi1,FMIN,ERRNUM] = SolveWithGLPK(H,phi0',M);

# Find Kantorovich min for ground state phi0 using 'linprog'
msg("Find gamma2 ...");
[gamma2,mu2,nu2,x,k,EK2,ES2,HM,phi2,FMIN] = SolveWithLinprog(H,phi0',M);

# Check results
msg("Check results ...");
checkGLPKSolution(gamma1,mu1,nu1,EK1,ERRNUM,FMIN);
checkGLPKSolution(gamma2,mu2,nu2,EK2,0,EK2);

# Report



# Plot
#plot(x,phi0)
#contour(x,k,gamma1)
#contourf(x,k,tensorProductMeasure(mu1,nu1),24)
# plot3(x,k,gamma1)

msg("=========================");
rep("ES0 .........",ES0);
rep("ES1 .........",ES1);
rep("ES2 .........",ES2);
rep("EK1 .........",EK1);
rep("EK2 .........",EK2);
msg("=========================");


