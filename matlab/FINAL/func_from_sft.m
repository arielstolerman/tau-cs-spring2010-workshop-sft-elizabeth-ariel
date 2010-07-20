function[res] = func_from_sft(L,coeffs,x,G)
res = 0;
size = length(G);
for ind=1:length(L)
    alpha = L(ind:ind,1:size);
    chi = 1;
    for j=1:size
        t = 2*pi./G(j)*alpha(j)*x(j);
        chi = chi*(cos(t)+i*sin(t));
    end
    res = res + coeffs(ind)*chi;
end