sft_setenv
isLogged = true;
G = 10^10;
tau = 200;
func = @(x,G)test(x,G);
numOfIterations = 1;
delta_t = 0.01;
fInfNorm = 28.41;
fEucNorm = 20;
deltaCoeff = 1;
maCoeff = 0.0001;
mbCoeff = 0.0001;
etaCoeff = 1;
[L,coeffs] = sft_dp_full(isLogged,G,tau,func,numOfIterations,delta_t,fInfNorm,fEucNorm,deltaCoeff,maCoeff,mbCoeff,etaCoeff)
