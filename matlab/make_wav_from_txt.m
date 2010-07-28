function[] = make_wav_from_txt(name)

% read from sft text output
f = fopen(strcat(name,'.txt'));
a = fscanf(f,'%g %g',[2,inf]);
fclose(f);
a = transpose(a);

% create normalized function in a
m = max(max(a));
for i=1:size(a,1)
   a(i,1) = a(i,1)./m; 
   a(i,2) = a(i,2)./m;
end

% write to wav files
wavwrite(a,44100,16,strcat(name,'.wav'));

