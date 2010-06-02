function[res]=test(x,G)
% testing a function over G -> C for a finite Abelian domain
% this is a simple example ment for G being {(1,N)} for some N

alpha = [230 1492 542 29];
coeff_alpha = [(10+10i) (10+10i) (0.01+0.01i) (0.07+0.09i)];

Gd = G(2);
res = 0;
for j=1:4
    term = 2*pi*alpha(j)*x/Gd;
    res = res + coeff_alpha(j)*(cos(term)+i*sin(term));
end
