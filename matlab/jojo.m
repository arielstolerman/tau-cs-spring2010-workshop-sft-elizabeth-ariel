function[res]=jojo(x,G)
% this is an automatically generated function from a FourierPolynomial java object.

k = length(G);
alpha = [1492 800 9221;230 500 28;29 29 9238;542 49 2918];
alpha = transpose(alpha);
coeff_alpha = [(19.0+120.0i) (101.0+130.0i) (0.07+0.09i) (0.01+0.01i) ];

res = 0;
size = size(alpha); alpha_size = size(1)*size(2);
for j=1:k:alpha_size;
vec = alpha(j:(k-1));	tmp = 1;
	for l=1:k;
		curr_alpha = vec(l);
		term = 2*pi*curr_alpha*x(l)/G(l);
		tmp = tmp * coeff_alpha(l)*(cos(term)+i*sin(term));
	end
res = res + tmp;
end
