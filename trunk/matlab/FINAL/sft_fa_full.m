% The SFT algorithm implementation for MATLAB over finite Abelian group domain
% Parameters:
% @param isLogged:		A flag to indicate weather to log the algorithm actions or not.
% @param G			The values (g1,N1),...,(gk,Nk) describing the Abelian group G where gj are the
%	 				corresponding generators for Nj. G should be a [2 X k] matrix, where each column is a vector transpose[gj Nj].
% @param tau	 		The threshold such that all tau-significant elements are returned. 
% @param func			The given function over G -> C whose tau-significant elements and thier Fourier coefficients are returned. Used for query access.
%					func should be a MATLAB function with parameters x (a vector in G) and G (a vector of pairs of g's and N's).
% @param numOfIterations	The number of SFT procedure iterations to run. Each iteration is ran with the difference function
%	 				of the given function and the output of the previous SFT iteration.
%	 				This is an optimization for the original SFT algorithm to enable catching significant coefficients
%	 				with greater precision.
% @param delta_t			The confidence parameter such that the algorithm succeeds with probability 1-delta.
% @param fInfNorm		The infinity norm of the function.
% @param fEuclideanNorm	The Euclidean norm of the function.
% @param deltaCoeff		A constant coefficient for the algorithm's calculation of delta.
% @param maCoeff		A constant coefficient for the algorithm's calculation of m_A.
% @param mbCoeff		A constant coefficient for the algorithm's calculation of m_B.
% @param etaCoeff		A constant coefficient for the algorithm's calculation of eta (appears in m_A and m_B calculation).
% Result:
% Returns a mapping of the elements in G and their tau-significant coefficients in the given function with delta-confidence.
% L - a vector of the tau-significant elements.
% coeffs - a vector of their corresponding coefficients (s.t. the coefficient of L(i) is coeffs(i)).

function[L,coeffs]=sft_fa_full(isLogged,G,tau,func,numOfIterations,delta_t,fInfNorm,fEuclideanNorm,deltaCoeff,maCoeff,mbCoeff,etaCoeff);

% set java path
import java.util.*
import java.io.File
import SFT.*
import SFT.SFTUtils.*

% calculate the corresponding direct product G and function
dim = size(G,2);
dpG = zeros(dim);
for i=1:dim
	dpG(i) = G(2,i);
end

% instead of the given 2nd argument G, which is to be given a direct product G, the function will always use the finite Abelian G
dpfunc = @(x,nonce)dp_func_res_from_fa_func(G,func,x);
[dpL,dpCoeffs] = sft_dp_full(isLogged,dpG,tau,dpfunc,numOfIterations,delta_t,fInfNorm,fEuclideanNorm,deltaCoeff,maCoeff,mbCoeff,etaCoeff);

% calculate the corresponding finite Abelian L
L = zeros(size(dpL,1),1);
for i=1:size(dpL,1)
	x = dpL(i);
	fax = 1;
	for j=1:dim
		gj = G(1,j);
		Nj = G(2,j);
		xj = x(j);
		fax = fax*mod(gj*xj,Nj);
	end
	L(i) = fax;
end

% coefficients remain the same
coeffs = dpCoeffs;
