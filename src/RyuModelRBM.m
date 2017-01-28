function [ERROR, J] = RyuModelRBM(UNKNOWNS, alpht, gamt, alphn, gamn, ds, vt, vn, tx, ty, nx, ny, x, y, xcm, ycm)

% 
% 
% The MIT License
% 
% Copyright (c) Eric Keaveny and André Brown
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.


Nx = size(x,2);

VRBMX = UNKNOWNS(1);
VRBMY = UNKNOWNS(2);
OMEG = UNKNOWNS(3);

ERROR = zeros(1, 3);
ft = 0.0*tx;
fn = 0.0*nx;
tor_int = 0.0*nx;

J = zeros(3,3);

for i = 1:Nx
    j = Nx + i;
    XmXcm = x(i) - xcm;
    YmYcm = y(i) - ycm;
    OmegCrossXx = -OMEG*YmYcm;
    OmegCrossXy = OMEG*XmXcm;
    
    OmegxXt = OmegCrossXx.*tx(i) + OmegCrossXy.*ty(i);
    OmegxXn = OmegCrossXx.*nx(i) + OmegCrossXy.*ny(i);
    
    vrbmt = VRBMX*tx(i) + VRBMY*ty(i);
    vrbmn = VRBMX*nx(i) + VRBMY*ny(i);
    
    vttot = vrbmt + OmegxXt + vt(i);
    vntot = vrbmn + OmegxXn + vn(i);
    
    ft(i) = alpht*vttot*abs(vttot).^(gamt - 1);
    fn(i) = alphn*vntot*abs(vntot).^(gamn - 1);
    tor_int(i) = XmXcm*(ft(i)*ty(i) + fn(i)*ny(i)) - YmYcm*(ft(i)*tx(i) + fn(i)*nx(i));
    
    dftdvt = alpht*gamt*abs(vttot).^(gamt - 1);
    dfndvn = alphn*gamn*abs(vntot).^(gamn - 1);
    
    dvtdvx = tx(i);
    dvtdvy = ty(i);
    dvtdw = XmXcm*ty(i) - YmYcm*tx(i);
    
    dvndvx = nx(i);
    dvndvy = ny(i);
    dvndw = XmXcm*ny(i) - YmYcm*nx(i);
    
    J(1,1) = J(1,1) + ds*dftdvt*dvtdvx*tx(i) + ds*dfndvn*dvndvx*nx(i);
    J(1,2) = J(1,2) + ds*dftdvt*dvtdvy*tx(i) + ds*dfndvn*dvndvy*nx(i);
    J(1,3) = J(1,3) + ds*dftdvt*dvtdw*tx(i) + ds*dfndvn*dvndw*nx(i);
    
    J(2,1) = J(2,1) + ds*dftdvt*dvtdvx*ty(i) + ds*dfndvn*dvndvx*ny(i);
    J(2,2) = J(2,2) + ds*dftdvt*dvtdvy*ty(i) + ds*dfndvn*dvndvy*ny(i);
    J(2,3) = J(2,3) + ds*dftdvt*dvtdw*ty(i) + ds*dfndvn*dvndw*ny(i);
    
    J(3,1) = J(3,1) + ds*dftdvt*dvtdvx*(XmXcm*ty(i) - YmYcm*tx(i)) + ds*dfndvn*dvndvx*(XmXcm*ny(i) - YmYcm*nx(i));
    J(3,2) = J(3,2) + ds*dftdvt*dvtdvy*(XmXcm*ty(i) - YmYcm*tx(i)) + ds*dfndvn*dvndvy*(XmXcm*ny(i) - YmYcm*nx(i));
    J(3,3) = J(3,3) + ds*dftdvt*dvtdw*(XmXcm*ty(i) - YmYcm*tx(i)) + ds*dfndvn*dvndw*(XmXcm*ny(i) - YmYcm*nx(i));
    
end

ERROR(1) = ds*sum(ft.*tx + fn.*nx);
ERROR(2) = ds*sum(ft.*ty + fn.*ny);
ERROR(3) = ds*sum(tor_int);
