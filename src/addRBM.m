function [VXmod, VYmod, X, Y] = addRBM(DX, DY, VX, VY, RBM, dt)

% ADDRBM adds predicted rigid body motion onto a series of skeletons and
% outputs the resulting velocities and positions.
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


Nt = size(VX,1);
Nseg = size(VX,2);

VXmod = zeros(Nt,Nseg);
VYmod = zeros(Nt,Nseg);

% add rigid body motion (xy speeds and rotation speed) to VX and VY
for ii=1:Nt
    VXmod(ii,:) = VX(ii,:) + RBM(ii,1) - RBM(ii,3).*DY(ii+1,:);
    VYmod(ii,:) = VY(ii,:) + RBM(ii,2) + RBM(ii,3).*DX(ii+1,:);
end

% use new speeds to calculate new X and Y coordinates
X = NaN(size(DX));
Y = NaN(size(DY));

% assign starting coordinates
X(1, :) = DX(1, :);
Y(1, :) = DY(1, :);

% use velocity to shift subsequent skeletons
for ii = 2:Nt
    X(ii, :) = X(ii-1, :) + VXmod(ii-1, :)*dt;
    Y(ii, :) = Y(ii-1, :) + VYmod(ii-1, :)*dt;
end