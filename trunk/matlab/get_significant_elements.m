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
Garr= javaArray('java.lang.long',length(G));
sft = SFT.SFT();
sets = sft.runMatlabSFTPart1Internal(Garr,delta_t,tau,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff);

tmp=sets(sets.length);
q= tmp(1).toArray;

query=javaObject('java.util.HashMap');

for i=1:q.length;
  xLong=q(i);
%  xLong=javaObject('java.lang.Long',x);
    x= zeros(1, xLong.length);
    for j=1:xLong.length;
        x(j)=xLong(j);
    end
  y=func(x,G);
  yComplex = Complex(real(y),imag(y));
  %yComplex=xmlFunc.getValue(x);
  query.put(xLong,yComplex);
end

L=sft.runMainSFTAlgorithmDividedPart2(Garr,tau,sets,query);
