% The SFT algorithm implementation for MATLAB over direct group-product domain(i.e. Z_N1 x ... x Z_Nk).
% Parameters:
% @param isLogged:		A flag to indicate weather to log the algorithm actions or not.
% @param G			The values N1,...,Nk describing the group G = Z_N1 x ... x Z_Nk.
% @param tau	 		The threshold such that all tau-significant elements are returned. 
% @param func			The given function over G -> C whose tau-significant elements and thier Fourier coefficients are returned. Used for query access.
%					func should be a MATLAB function with parameters x (a vector in G) and G (a vector of N's).
% @param m_A			The size of the group A (for constructing the group Q).
% @param m_B			The size of the groups Btl (t=1,...,k, l=1,...,log2(Nt), for constructing the group Q).
% @param numOfIterations	The number of SFT procedure iterations to run. Each iteration is ran with the difference function
%	 				of the given function and the output of the previous SFT iteration.
%	 				This is an optimization for the original SFT algorithm to enable catching significant coefficients
%	 				with greater precision.
% Result:
% Returns a mapping of the elements in G and their tau-significant coefficients in the given function with confidence set by the values m_A and m_B.
% L - a vector of the tau-significant elements.
% coeffs - a vector of their corresponding coefficients (s.t. the coefficient of L(i) is coeffs(i)).

function[L,coeffs]=sft_dp(isLogged,G,tau,func,m_A,m_B,numOfIterations);

if nargin == 6;
	numOfIterations = 1;
end;

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
rep = sft.runMatlabSFTPart1Internal(isLogged_java,G_java,tau,m_A,m_B);

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
jrandSet = jres.getRandSet;
randSetSize = jrandSet.size;
L = zeros(jkeys.length,size);	% for holding the significant elements
coeffs = zeros(jkeys.length,1);	% for holding their coefficients
for ind=1:jkeys.length;
	xLong=jkeys(ind);
  	for j=1:size;
    	L(ind,j)=xLong(j).longValue;
  	end
	x = L(ind,1:size);
	val = jvalues(ind);
	% calculate the coefficient over the random set of elements jrandSet
	coeffTmp = complex(0,0);
	for k=1:randSetSize
		yLong=jrandSet(k);
		y = 0;
		for j=1:size;
			y(j) = yLong(j).longValue;
		end
		term = 2*pi./G;
		chi = 1;
		for j=1:size
			tmpTerm = term(j)*y(j)*x(j);
			chi = chi*complex(cos(tmpTerm),sin(tmpTerm));
		end
		chi = conj(chi);
		coeffTmp = coeffTmp + func(x,G)*chi./randSetSize;
	end
  	coeffs(ind) = coeffTmp;
end