% The SFT algorithm implementation for MATLAB
% Parameters:
% - G: 				Vector of Integers representing Z_N1 x ... x Z_Nk
% - delta_t: 		The confidence parameter such that the algorithm succeeds with probability 1-delta.
% - tau:			The threshold such that all tau-significant elements are returned.
% - func:			the function for querying, over G -> C
% - fInfNorm:		the infinity norm of the function
% - fEuclideanNorm:	the Euclidean norm of the function
% - deltaCoeff:		a constant coefficient for calculating delta
% - randSetsCoeff:	a constant coefficient for calculating the random subsets when creating Q, a set of elements for querying
% - isLogged:		a flag to indicate weather to log the algorithm actions or not
% TODO: add description

function[L]=get_significant_elements(G,delta_t,tau,func,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff,isLogged);

javaaddpath('/specific/a/home/cc/students/cs/arielst1/sft/sft_lib.jar')
import java.util.*
import java.io.File
import SFT.*
import SFT.SFTUtils.*

% comply all parameters to java before calling the first part of the algorithm
G_java = javaArray('java.lang.Long',length(G));
for i=1:length(G);
    G_java(i) = javaObject('java.lang.Long',G(i));
end
isLogged_java = javaObject('java.lang.Boolean',isLogged);	

% create a static object in order to allow calls to the SFT algorithm methods
sft = SFT.SFT();

% ===================
% =		PART 1		=
% ===================
% call part 1 - create a set Q of elements (vectors) in G to be queried
rep = sft.runMatlabSFTPart1Internal(G_java,delta_t,tau,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff,isLogged_java);

% create query
q = rep.getQ;
query=javaObject('java.util.HashMap');

% calculate function values on Q's elements
utils = SFT.SFTUtils(); 
for i=1:q.length;
  xLong=q(i);
  x=zeros(1, length(xLong));
  for j=1:length(xLong);
      x(j)=xLong(j).longValue;
  end
  y=func(x,G);
  yComplex = Complex(real(y),imag(y));
  query.put(utils.vectorToString(xLong),yComplex);
end
% update temporary repository
rep.setQuery(query);

% ===================
% =		PART 2		=
% ===================
% call part 2 - create an array of Long[] from which will create a final matlab matrix where
% each row is an element (vector) in L, the set of significant elements
L_java=sft.runMatlabSFTPart2Internal(G_java,tau,rep);
L = zeros(L_java.length,length(G));
for i=1:L_java.length;
	xLong=L_java(i);
	x=zeros(1, length(xLong));
  	for j=1:length(xLong);
    	x(j)=xLong(j).longValue;
  	end
  	L(i)=x;							%% WILL NOT WORK
end

