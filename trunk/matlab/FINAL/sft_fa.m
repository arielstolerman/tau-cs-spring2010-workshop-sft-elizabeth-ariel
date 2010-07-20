% The SFT algorithm implementation for MATLAB over finite Abelian group domain
% Parameters:
% @param isLogged:		A flag to indicate weather to log the algorithm actions or not.
% @param G			The values (g1,N1),...,(gk,Nk) describing the Abelian group G where gj are the
%	 				corresponding generators for Nj. G should be a [2 X k] matrix, where each column is a vector transpose[gj Nj].
% @param tau	 		The threshold such that all tau-significant elements are returned. 
% @param func			The given function over G -> C whose tau-significant elements and thier Fourier coefficients are returned. Used for query access.
%					func should be a MATLAB function with parameters x (a vector in G) and G (a vector of pairs of g's and N's).
% @param m_A			The size of the group A (for constructing the group Q).
% @param m_B			The size of the groups Btl (t=1,...,k, l=1,...,log2(Nt), for constructing the group Q).
% @param numOfIterations	The number of SFT procedure iterations to run. Each iteration is ran with the difference function
%	 				of the given function and the output of the previous SFT iteration.
%	 				This is an optimization for the original SFT algorithm to enable catching significant coefficients
%	 				with greater precision.
% Result:
% Returns a mapping of the elements in G and their tau-significant coefficients in the given function with delta-confidence.
% L - a vector of the tau-significant elements.
% coeffs - a vector of their corresponding coefficients (s.t. the coefficient of L(i) is coeffs(i)).

function[L,coeffs]=sft_fa(isLogged,G,tau,func,m_A,m_B,numOfIterations);

% set java path
import java.util.*
import java.io.File
import SFT.*
import SFT.SFTUtils.*

% fit all parameters to java before calling the first part of the algorithm
G_java = javaArray('java.lang.Long[]',size(G,2));
for i=1:size(G,2);
	tmp = javaArray('java.lang.Long',2);
	tmp(1) = javaObject('java.lang.Long',G(1,i));
	tmp(2) = javaObject('java.lang.Long',G(2,i));
	G_java(i) = tmp;
end
isLogged_java = javaObject('java.lang.Boolean',isLogged);	

% create a static object in order to allow calls to the SFT algorithm methods
sft = SFT.SFT();

% ======
% PART 1
% ======
% call part 1 - create a set Q of elements (vectors) in G to be queried
rep = sft.runMatlabSFTPart1Internal(isLogged_java,G_java,tau,m_A,m_B);

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
  query.put(xLong.toString,yComplex);
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
L = zeros(jkeys.length,1);		% for holding the significant elements
coeffs = zeros(jkeys.length,1);	% for holding their coefficients
for ind=1:jkeys.length;
	xLong=jkeys(ind);
    L(ind)=xLong.longValue;
	val = jvalues(ind);
  	coeffs(ind) = complex(val(1).doubleValue,val(2).doubleValue);
end