% defines a function over the matrix m
% where m is an output of the Matlab function WAVREAD
% in order to use a wav file as an input function for the SFT scripts, do the following:
%	m = wavread('<WAV file path>');
%	func = @(x,G)wav_func_script(x,m);
% (note that G has no affect in the definition above);

function[res] = wav_func_script(x,m);
tmp = m(x+1,1:2);
res = complex(tmp(1),tmp(2));