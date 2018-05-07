function str = gen_ker_estrin_C_code(w,d,be,o)
% GEN_KER_ESTRIN_C_CODE  Write C code strings for Estrin eval of ES kernel
%
% str = gen_ker_estrin_C_code(w,d,be,o)
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
% but that might be a bit faster; not sure. Stuck to simple Estrin for now,
% not exploiting that there are w calls to different poly's w/ *same* z arg.

% Barnett 4/23/18.
% Ludvig af Klinteberg, 7 May 2018

if nargin==0, test_gen_ker_estrin_C_code; return; end
if nargin<4, o=[]; end

C = ker_ppval_coeff_mat(w,d,be,o);
str = cell(d+1,1);

% Coeff arrays
for n=1:d % loop over poly coeffs
    s = sprintf('FLT c%d[] = {%.16E',n-1, C(n,1));
    for i=2:w % loop over segments
        s = sprintf('%s, %.16E', s, C(n,i));      
    end
    str{n} = [s sprintf('};\n')];
end

s = '';
nl = sprintf('\n');
tab = '  ';
tab2 = [tab tab];
tab3 = [tab2 tab];

% Powers of z that are needed
n = 2;
s = [s sprintf('FLT z2 = z*z;\n') tab2];
while n*2<d
    s = [s sprintf('FLT z%d = z%d*z%d;\n', 2*n, n, n) tab2];
    n = 2*n;
end

% Recursively form code for Estrin's scheme
ops = estrin_line(d, 0, d);

s = [s sprintf('for (int i=0; i<%d; i++)\n',w) tab3];
s = [s 'ker[i] = ' nl tab3];
s = [s ops ';' nl];

str{d+1} = s;
end

function [s, count] = estrin_line(n, count, max)
% Recursively form powers up to z^n using Estrin's scheme
%
% n = order of base polynomial
% count = number of next poly coeff to use
% max = largest poly coeff (i.e. final order)
if n==1
    if count+1<max
        % even poly order
        s = sprintf('c%d[i] + c%d[i]*z', count, count+1);
        count = count+2;
    else 
        % odd
        s = sprintf('c%d[i]', count);
        count = count+1;
    end        
else
    m = 1;
    s = '';
    while m<n && count<max
        [p, count] = estrin_line(m, count, max);
        if m==1
            s = sprintf('%s%s',s, p);                        
        else
            s = sprintf('%s + z%d*(%s)',s, m, p);            
        end
        m = 2*m;        
    end
end
end

%%%%%%%%
function test_gen_ker_estrin_C_code    % writes C code to file but doesn't test.
w=13; d=16;           % pick a single kernel width and degree to write code for
                      %w=7; d=11;
                      %w=2; d=5;
beta=2.3*w;
str = gen_ker_estrin_C_code(w,d,beta);
% str{:}
fnam = sprintf('ker_estrin_w%d.c',w);
fid = fopen(fnam,'w');
for i=1:numel(str); fwrite(fid,str{i}); end
fclose(fid);
system(['more ' fnam])
end