% Based on function polywind from NR, but vectorized
% One change: input vectors will already by closed, so we just run the full
% loop rather than starting at the "last" vertex.

function wind = polywindFlag(vt, pt, flag)

if nargin < 3 || flag == 0
    pt0 = pt(1);
    pt1 = pt(2);
    
    d0 = vt(2:end, 1);
    d1 = vt(2:end, 2);
    p0 = vt(1:end-1, 1);
    p1 = vt(1:end-1, 2);
    
    windp = (p1 <= pt1) & (d1 >  pt1) & (((p0 - pt0).*(d1-pt1)-(p1-pt1).*(d0-pt0)) > 0);
    windn = (p1 > pt1) & (d1 <= pt1) & (((p0 - pt0).*(d1-pt1)-(p1-pt1).*(d0-pt0)) < 0);
    wind = sum(windp) - sum(windn);
elseif flag == 1 % all points known a priori to be inside occulter.
    wind = 1;
end
