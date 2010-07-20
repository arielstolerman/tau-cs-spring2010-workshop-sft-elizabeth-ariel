% given an output of the SFT algorithm over a WAV file (L and coeffs),
% the domain N and a desired name, this script creates a WAV with that name

function[] = make_wav_from_sft(L,coeffs,N,name)

% calculate all values for 0,...,N-1
t = zeros(N,2);
for i=1:N
    if (mod(i,100) == 0)
        fprintf(1,'done %d out of %d calculations\n',i,N)
    end
    tmp = func_from_sft(L,coeffs,i,N);
    t(i,1) = real(tmp);
    t(i,2) = imag(tmp);
end
% normalize the matrix
m1 = max(max(t));
m2 = abs(min(min(t)));
m = max(m1,m2);
t = t./m;
% create WAV file
wavwrite(t,44100,16,strcat(name,'.wav'));