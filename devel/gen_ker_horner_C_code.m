function str = gen_ker_horner_C_code(w,d,be,o)
% GEN_KER_HORNER_C_CODE  Write C code strings for Horner eval of ES kernel
%
% str = gen_ker_horner_C_code(w,d,be,o)
%
% Inputs:
%  w = integer kernel width in grid points, eg 10
%  d = poly degree to keep, eg 13
%  beta = kernel parameter, around 2.3*w
%  opts - optional struct (unused; could switch to cosh kernel variant, etc..)
%
% Outputs:
%  str = length-w cell array of C code strings to eval each segment of kernel
%
% Also see: KER_PPVAL_COEFF_MAT, FIG_SPEED_KER_PPVAL (which tests acc too)
%
% Note: # flops is same as filling a col vec of [1;z;z^2;..] & doing small BLAS2
% but that might be a bit faster; not sure. Stuck to simple Horner for now,
% not exploiting that there are w calls to different poly's w/ *same* z arg.

% Barnett 4/23/18.
if nargin==0, test_gen_ker_horner_C_code; return; end
if nargin<4, o=[]; end

C = ker_ppval_coeff_mat(w,d,be,o);
str = cell(d+1,1);
for n=1:d % loop over poly coeffs
  s = sprintf('FLT c%d[] = {%.16E',n-1, C(n,1));
  for i=2:w % loop over segments
    s = sprintf('%s, %.16E', s, C(n,i));      
  end
  str{n} = [s sprintf('};\n')];
end

s = sprintf('for (int i=0; i<%d; i++) ker[i] = ',w);
for n=1:d-1
  s = [s sprintf('c%d[i] + z*(',n-1)];   % (n-1)th coeff for i'th segment
end
s = [s sprintf('c%d[i]',d-1)];
for n=1:d-1, s = [s sprintf(')')]; end  % close all parens
s = [s sprintf(';\n')];          % terminate the C line, CR
str{d+1} = s;

%%%%%%%%
function test_gen_ker_horner_C_code    % writes C code to file but doesn't test.
w=13; d=16;           % pick a single kernel width and degree to write code for
%w=7; d=11;
%w=2; d=5;
beta=2.3*w;
str = gen_ker_horner_C_code(w,d,beta);
% str{:}
fnam = sprintf('ker_horner_w%d.c',w);
fid = fopen(fnam,'w');
for i=1:numel(str); fwrite(fid,str{i}); end
fclose(fid);
system(['more ' fnam])
