function RBM = posture2RBM_nonlinear(X, Y, XCM, YCM, TX, TY, VT, VN, NX, NY, ...
    ds, alpht, alphn, gamt, gamn)

% POSTURE2RBM_nonlinear converts skeleton information and drag coefficients
% ratio into a matrix of rigid body motions.
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


% initialise
RBM = zeros(size(X, 1) - 1, 3);
FrameRBM = ones(1,3);

% set the options for fsolve
fsolveopts = optimset('Display','off','Jacobian','on', ...
    'MaxFunEvals',1000,'MaxIter',1000, 'TolFun', 1e-12);

% loop through frames
for ii = 1:size(X, 1) - 1    
    % run the model
    ryumod_funhan = @(UNKNOWNS)RyuModelRBM(UNKNOWNS, ...
        alpht, gamt, alphn, gamn, ds, VT(ii,:), VN(ii,:), ...
        TX(ii,:), TY(ii,:), NX(ii,:), ...
        NY(ii,:), X(ii,:), Y(ii,:), XCM(ii), YCM(ii));
    FrameRBM = fsolve(ryumod_funhan, FrameRBM, fsolveopts);  % first few frames may be innaccurate because starting point is just ones
    RBM(ii,:) = FrameRBM(1:end);
end