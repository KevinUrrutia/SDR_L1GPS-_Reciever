function [I, Q] = IQ_parsing(fname)
%inputs: file name-> string containing file path
%outputs sampled real and imaginary values from the sampled data

fid = fopen(fname, 'r', 'l');
fseek(fid, 0, 'eof');
N = (ftell(fid) / 4) / 2; %#bytes* 2samples/4bytes * sample/2
fseek(fid, 0, "bof");

dat = fread(fid, [2, N], 'int16');
fclose(fid);
I = dat(1, :);
Q = dat(2, :);
end