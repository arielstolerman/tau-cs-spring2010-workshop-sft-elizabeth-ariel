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

function[resL,resCoeffs]=sft_dp(isLogged,G,tau,func,m_A,m_B,numOfIterations);

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

% ======
% PART 2
% ======
%  run numOfIterations iterations. In each iteration i:
% - calculate Q's values for function f_i
% - run Java part 2 of the SFT to get the list of significant elements and a random set to estimate the coefficients for those elements over it
% - calculate the coefficients for the given elements
% - calculate the next iteration function f_(i+1) = f - g_i, where f is the original input function and g_i is the SFT estimation for the i iteration
% The output: the union of all elements received from all iterations (and their coefficients).

% initialize first iteration's function and result L,coeffs
currFunc = func;
dim = length(G);
resL = zeros(1,dim);	% for holding the significant elements
resCoeffs = zeros(1,1);	% for holding their coefficients

% run iterations
for iter=1:numOfIterations
	if (isLogged)
		fprintf(1,'>>>>>>> MATLAB: starting iteration %d out of %d\n',iter,numOfIterations)
	end
	
	% create query
	q = rep.getQ;
	query=javaObject('java.util.HashMap');

	% calculate current iteration's function values on Q's elements
	utils = SFT.SFTUtils(); 
	for i=1:q.length;
	  xLong=q(i);
	  x=zeros(1, dim);
	  for j=1:dim;
		  x(j)=xLong(j).longValue;
	  end
	  y=currFunc(x,G);
	  yComplex = Complex(real(y),imag(y));
	  query.put(utils.vectorToString(xLong),yComplex);
	end
	% update temporary repository
	rep.setQuery(query);

	% call part 2 - create an array of Long[] from which will create a final matlab matrix where
	% each row is an element (vector) in L, the set of significant elements
	jres=sft.runMatlabSFTPart2Internal(G_java,tau,rep,1);
	
	% if no new elements caught for this iteration (and it is not the first iteration), break
	if (iter > 1 && jres.getKeys.length == 0)
		if (isLogged)
			fprintf(1,'>>>>>>> MATLAB: iteration %d finished with 0 new elements, breaking calculation.\n',iter)
		end
		break;
	end

	jkeys = jres.getKeys;
	jrandSet = jres.getRandSet;
	randSetSize = jrandSet.length;
	L = zeros(jkeys.length,dim);	% for holding the significant elements
	coeffs = zeros(jkeys.length,1);	% for holding their coefficients
	if (isLogged)
		fprintf(1,'>>>>>>> MATLAB: starting coefficients calculation in iteration %d\n',iter);
	end
	for ind=1:jkeys.length;
		if (isLogged && mod(ind,100) == 0)
			fprintf(1,'Done %d coefficients calculations...\n',ind);
		end
		
		xLong=jkeys(ind);
		for j=1:dim;
			L(ind,j)=xLong(j).longValue;
		end
		x = L(ind,1:dim);

		% calculate the coefficient over the random set of elements jrandSet
		coeffTmp = complex(0,0);
		for k=1:randSetSize
			yLong=jrandSet(k);
			y = zeros(dim);
			for j=1:dim;
				y(j) = yLong(j).longValue;
			end
			chi = 1;
			for j=1:dim
				term = 2*pi*(y(j)./G(j))*x(j);
				chi = chi*complex(cos(term),sin(term));
			end
			chi = conj(chi);
			coeffTmp = coeffTmp + (func(y,G)*chi)./randSetSize;
		end
		coeffs(ind) = coeffTmp;
	end
	
	% initialize result L,coeffs if this is the first iteration
	if (iter == 1)
		resL = L;
		resCoeffs = coeffs;
	% otherwise only update resL and resCoeffs
	else
		tmpResL = resL;
		tmpResCoeffs = resCoeffs;
		sizeOfRes = size(resL,1);
		sizeOfCurr = size(L,1);
		ind = sizeOfRes + 1; % index to start adding new elements from
		for i=1:sizeOfCurr
			isContained = false;
			for j=1:sizeOfRes
				if (min(L(i) == resL(j)))
					isContained = true;
					break;
				end
			end
			if (~isContained)
				% add element and coeff to final result
				tmpResL(ind,1:dim) = L(i,:);
				tmpResCoeffs(ind) = coeffs(i);
				ind = ind + 1;
			end
		end
		resL = tmpResL;
		resCoeffs = tmpResCoeffs;
	end
	
	% update function for next iteration
	currFunc = @(x,G)(func(x,G) - func_from_sft(resL,resCoeffs,x,G));
	
	if (isLogged)
		fprintf(1,'>>>>>>> MATLAB: finished iteration %d out of %d.\n',iter,numOfIterations)
	end
end % iteration