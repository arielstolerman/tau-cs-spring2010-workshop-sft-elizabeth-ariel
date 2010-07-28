% given a finiteAbelian domain representation faG, a function over this domain
% and an element x in the corresponding direct product domain, returns the function's value over the corresponding
% finite Abelian element
function[res] = dp_func_res_from_fa_func(faG,faFunc,x)
% calculate the finite Abelian element
dim = size(faG,2);
fax = 1;
for i=1:dim
	gi = faG(1,i);
	Ni = faG(2,i);
	xi = x(i);
	fax = fax*mod(gi*xi,Ni);
end
res = faFunc(fax,faG);
