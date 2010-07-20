% The SFT algorithm implementation for MATLAB over direct group-product domain(i.e. Z_N1 x ... x Z_Nk).
% Parameters:
% @param isLogged:		A flag to indicate weather to log the algorithm actions or not.
% @param G			The values N1,...,Nk describing the group G = Z_N1 x ... x Z_Nk.
% @param tau	 		The threshold such that all tau-significant elements are returned. 
% @param func			The given function over G -> C whose tau-significant elements and thier Fourier coefficients are returned. Used for query access.
%					func should be a MATLAB function with parameters x (a vector in G) and G (a vector of N's).
% @param numOfIterations	The number of SFT procedure iterations to run. Each iteration is ran with the difference function
%	 				of the given function and the output of the previous SFT iteration.
%	 				This is an optimization for the original SFT algorithm to enable catching significant coefficients
%	 				with greater precision.
% @param delta_t			The confidence parameter such that the algorithm succeeds with probability 1-delta.
% @param fInfNorm		The infinity norm of the function.
% @param fEuclideanNorm	The Euclidean norm of the function.
% @param deltaCoeff		A constant coefficient for the algorithm's calculation of delta.
% @param randSetsCoeff	A constant coefficient for the algorithm's calculation of the random sets A, Btl (t=1,...,k, l=1,...,log2(Nt)).
% Result:
% Returns a mapping of the elements in G and their tau-significant coefficients in the given function with delta-confidence.
% L - a vector of the tau-significant elements.
% coeffs - a vector of their corresponding coefficients (s.t. the coefficient of L(i) is coeffs(i)).

function[L,coeffs]=sft_dp_full(isLogged,G,tau,func,numOfIterations,delta_t,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff);

% set java path
import java.util.*
import java.io.File
import SFT.*
import SFT.SFTUtils.*

% fit all parameters to java before calling the first part of the algorithm
G_java = javaArray('java.lang.Long',length(G));
for i=1:length(G);
    G_java(i) = javaObject('java.lang.Long',G(i));
end
isLogged_java = javaObject('java.lang.Boolean',isLogged);	

% create a static object in order to allow calls to the SFT algorithm methods
sft = SFT.SFT();

% ======
% PART 1
% ======
% call part 1 - create a set Q of elements (vectors) in G to be queried
rep = sft.runMatlabSFTPart1Internal(isLogged_java,G_java,tau,delta_t,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff);

% create query
q = rep.getQ;
query=javaObject('java.util.HashMap');

% calculate function values on Q's elements
utils = SFT.SFTUtils(); 
size = length(G);
for i=1:q.length;
  xLong=q(i);
  x=zeros(1, size);
  for j=1:size;
      x(j)=xLong(j).longValue;
  end
  y=func(x,G);
  yComplex = Complex(real(y),imag(y));
  query.put(utils.vectorToString(xLong),yComplex);
end
% update temporary repository
rep.setQuery(query);

% ======
% PART 2
% ======
% call part 2 - create an array of Long[] from which will create a final matlab matrix where
% each row is an element (vector) in L, the set of significant elements
jres=sft.runMatlabSFTPart2Internal(G_java,tau,rep,numOfIterations);

jkeys = jres.getKeys;
jvalues = jres.getValues;
L = zeros(jkeys.length,size);	% for holding the significant elements
coeffs = zeros(jkeys.length,1);	% for holding their coefficients
for ind=1:jkeys.length;
	xLong=jkeys(ind);
  	for j=1:size;
    	L(ind,j)=xLong(j).longValue;
  	end
	val = jvalues(ind);
  	coeffs(ind) = complex(val(1).doubleValue,val(2).doubleValue);
end