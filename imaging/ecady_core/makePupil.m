% makePupil
% Eric Cady, last rev. 7/18/08
%
% Creates an open pupil with a secondary.

function pupil = makePupil(Nx, Ny, OR, IR, deltax, deltay)
dx = 1/Nx;
dy = 1/Ny;
x = -1/2+dx/2:dx:1/2-dx/2;
y = -1/2+dy/2:dy:1/2-dy/2;
[Y, X] = meshgrid(y-deltay/Ny, x-deltax/Nx);
pupil = (((X.*X + Y.*Y) < (OR/2)^2) & ((X.*X + Y.*Y) >= (IR/2)^2)) + 0;