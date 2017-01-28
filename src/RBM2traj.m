function [XCM, YCM] = RBM2traj(RBM, dt)

% RBM2TRAJ converts a predicted rigid body motion matrix to a centre of
% mass trajectory.
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


Ntcm = size(RBM,1);

%Center of mass position as a function of time
XCM = zeros(Ntcm,1);
YCM = zeros(Ntcm,1);

for ii = 2:Ntcm
    itm1 = ii - 1;
    XCM(ii) = XCM(itm1) + dt*RBM(itm1,1);
    YCM(ii) = YCM(itm1) + dt*RBM(itm1,2);
end