% The SFT algorithm implementation for MATLAB
% Parameters:
% N - the value representing Z_N
% func - the function for query access by the algorithm
% TODO: add description

function[L]=get_significant_elements(G,delta_t,tau,func,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff);

Garr= javaArray('java.lang.Long',length(G));
for i=1:length(G);
    Garr(i) = javaObject('java.lang.Long',G(i));
end
sft = SFT.SFT();
sets = sft.runMatlabSFTPart1Internal(Garr,delta_t,tau,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff);

tmp=sets(sets.length);
q= tmp(1).toArray;

% NOTE: the problem here is that q(i) is not of java type long[], so it
% doesn't work. should do something about it... maybe change the returned
% type to use Long[] and not long[]

query=javaObject('java.util.HashMap');

for i=1:q.length;
  xLong=q(i);
  x= zeros(1, length(xLong));
  for j=1:length(xLong);
      x(j)=xLong(j);
  end
  y=func(x,G);
  yComplex = Complex(real(y),imag(y));
  query.put(xLong,yComplex);
end

L=sft.runMatlabSFTPart2Internal(Garr,tau,sets,query);
