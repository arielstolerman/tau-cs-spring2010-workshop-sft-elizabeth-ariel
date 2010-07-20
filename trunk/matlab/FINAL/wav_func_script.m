% defines a function over the matrix m
function[res] = wav_func_script(x,m);
tmp = m(x+1,1:2);
res = tmp(1) + i*tmp(2);
