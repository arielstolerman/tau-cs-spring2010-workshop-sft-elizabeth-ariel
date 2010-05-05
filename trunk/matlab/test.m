function[res]=test(x,N)
% testing a function over Z_N -> C
alpha = [2031200 1492 542 29];
coeff_alpha = [(10+10i) (10+10i) (0.01+0.01i) (0.07+0.09i)];

res = 0;
for j=1:4
    term = 2*pi*alpha(j)*x/N;
    res = res + coeff_alpha(j)*(cos(term)+i*sin(term));
end
