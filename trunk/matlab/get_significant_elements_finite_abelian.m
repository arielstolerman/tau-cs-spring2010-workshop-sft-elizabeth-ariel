% The SFT algorithm implementation for MATLAB
% Parameters:
% - isLogged:		a flag to indicate weather to log the algorithm actions or not
% - G: 				Matrix of Integers representing an Abelian group by orders and their generators (g1,N1),...,(gk,Nk)
%					should be a kx2 matrix s.t. each column c_j holds c_j(1) = g_j, c_j(2) = N_j 			
% - delta_t: 		The confidence parameter such that the algorithm succeeds with probability 1-delta.
% - tau:			The threshold such that all tau-significant elements are returned.
% - func:			the function for querying, over G -> C
% - fInfNorm:		the infinity norm of the function
% - fEuclideanNorm:	the Euclidean norm of the function
% OPTIONAL:
% - deltaCoeff:		a constant coefficient for calculating delta
% - randSetsCoeff:	a constant coefficient for calculating the random subsets when creating Q, a set of elements for querying
% Returns a list of elements in G whose Fourier coefficients are tau-significant.

function[L]=get_significant_elements_finite_abelian(isLogged,G,delta_t,tau,func,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff);

%javaaddpath('/specific/a/home/cc/students/cs/arielst1/sft/sft_lib.jar')
import java.util.*
import java.io.File
import SFT.*
import SFT.SFTUtils.*

% fit all parameters to java before calling the first part of the algorithm
G_java = javaArray('java.lang.Long[]',length(G));
for i=1:length(G);
	tmp = javaArray('java.lang.Long',2);
	tmp(1) = javaObject('java.lang.Long',G(1,i));
	tmp(2) = javaObject('java.lang.Long',G(2,i));
	G_java(i) = tmp;
end
isLogged_java = javaObject('java.lang.Boolean',isLogged);	

% create a static object in order to allow calls to the SFT algorithm methods
sft = SFT.SFT();

% if this is the "short" call (without the constants deltaCoeff and randSetsCoeff), get them
if nargin == 7;
	deltaCoeff = sft.getDeltaCoeff();
	randSetsCoeff = sft.getRandSetsCoeff();
end;

% ===================
% =		PART 1		=
% ===================
% call part 1 - create a set Q of elements in G to be queried
rep = sft.runMatlabSFTPart1Internal(G_java,delta_t,tau,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff,isLogged_java);

% create query
q = rep.getQ;
query=javaObject('java.util.HashMap');

% calculate function values on Q's elements
utils = SFT.SFTUtils(); 
for i=1:q.length;
  xLong=q(i);
  x=xLong.longValue;
  y=func(x,G);
  yComplex = Complex(real(y),imag(y));
  query.put(utils.vectorToString(xLong),yComplex);
end
% update temporary repository
rep.setQuery(query);

% ===================
% =		PART 2		=
% ===================
% call part 2 - create an array of Long from which will create a final matlab matrix where
% each row is an element (vector) in L, the set of significant elements
L_java=sft.runMatlabSFTPart2Internal(G_java,tau,rep);
for i=1:L_java.length;
	xLong=L_java(i);
  	L(i)=xLong.longValue;
end

