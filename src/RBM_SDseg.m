function SDvel = RBM_SDseg(UX, UY, TX, TY, DX, DY, VX, VY, L, I, ds, alpha)

% RBM_SD estimates the rigid body motion from skeleton data and returns the
% error between the estimate and the actual motion.
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


% first make the prediction
RBM = posture2RBM(TX, TY, DX, DY, VX, VY, L, I, ds, alpha);

% Compute the velocity of each segment using the RBM predicted by the model

Nt = size(VX,1);
Nseg = size(VX,2);

UXmod = zeros(Nt,Nseg);
UYmod = zeros(Nt,Nseg);

for ii=1:Nt
    UXmod(ii,:) = VX(ii,:) + RBM(ii,1) - RBM(ii,3).*DY(ii+1,:);
    UYmod(ii,:) = VY(ii,:) + RBM(ii,2) + RBM(ii,3).*DX(ii+1,:);
end

UERROR = sum((UX - UXmod).^2 + (UY - UYmod).^2, 2);

% estimate the error
SDvel = mean(UERROR);

% disp([alpha, SDvel*1e5])
