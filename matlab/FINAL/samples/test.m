function[res]=test(x,G)
% testing a function over G -> C for a direct product of finite groups domain
alpha = [230 1492 542 29];
coeff_alpha = [(10+10i) (10+10i) (0.01+0.01i) (0.07+0.09i)];

res = 0;
for j=1:4
    term = 2*pi*alpha(j)*x/G;
    res = res + coeff_alpha(j)*(cos(term)+i*sin(term));
end
