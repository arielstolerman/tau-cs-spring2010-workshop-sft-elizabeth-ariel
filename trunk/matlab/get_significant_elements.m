% The SFT algorithm implementation for MATLAB
% Parameters:
% N - the value representing Z_N
% func - the function for query access by the algorithm
% TODO: add description

function[L]=get_significant_elements(G,delta_t,tau,func,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff);
import SFT.*;
%import Function.*;
%import java.io.File;
%xmlFile = File('test.xml');
%xmlFunc = XMLFourierPolynomial(xmlFile,N);

sft = SFT.SFT();
sets = sft.runMainSFTAlgorithmDividedPart1(N,delta_t,tau,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff);

q=sets(sets.length).toArray;
query=javaObject('java.util.HashMap');

for i=1:q.length;
  x=q(i);
  xLong=javaObject('java.lang.Long',x);
  y=func(x,N);
  yComplex = Complex(real(y),imag(y));
  %yComplex=xmlFunc.getValue(x);
  query.put(xLong,yComplex);
end

L=sft.runMainSFTAlgorithmDividedPart2(N,tau,sets,query);
