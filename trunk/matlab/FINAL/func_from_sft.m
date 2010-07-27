%  given the SFT output L, coeffs for an input function over G -> C
% this script calculates the OUTPUT function's value in x

function[res] = func_from_sft(L,coeffs,x,G)

res = 0;
dim = length(G);
for ind=1:length(L)
    alpha = L(ind:ind,1:dim);
    chi = 1;
    for j=1:dim
        t = 2*pi./G(j)*alpha(j)*x(j);
        chi = chi*complex(cos(t),sin(t));
    end
    res = res + coeffs(ind)*chi;
end