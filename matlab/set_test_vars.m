% for direct procdut and finite abelian domains
isLogged=false;
G_direct_prod = 10^10;
G_finite_abelian = [1; 10^10];
delta_t=0.01;
tau=200;
func_direct_prod=@(x,G)example_direct_prod(x,G);
func_finite_abelian=@(x,G)example_finite_abelian(x,G);
fInfNorm=28.41;
fEuclideanNorm=20;
deltaCoeff=1;
randSetsCoeff=0.0001;
