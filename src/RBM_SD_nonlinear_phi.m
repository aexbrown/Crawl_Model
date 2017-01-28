function SDvel = RBM_SD_nonlinear_phi(UXCM, UYCM, X, Y, XCM, YCM, TX, TY, ...
    VT, VN, NX, NY, ds, alpht, alphn, gamt, gamn, phi)

% RBM_SD estimates the rigid body motion from skeleton data and returns the
% error between the estimate and the actual motion.  NB it only takes into
% account the translational velocity prediction, not the angular velocity
% prediction
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


% scale parameters by psi
[alpht, alphn, gamt, gamn] = phiTransform(alpht, alphn, gamt, gamn, phi);

% first make the prediction
RBM = posture2RBM_nonlinear(X, Y, XCM, YCM, TX, TY, VT, VN, NX, NY, ...
    ds, alpht, alphn, gamt, gamn);

% estimate the error
SDvel = mean((RBM(:, 1) - UXCM).^2 + (RBM(:, 2) - UYCM).^2);

% disp([phi, SDvel*1e5])