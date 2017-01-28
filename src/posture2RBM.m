function RBM = posture2RBM(TX, TY, DX, DY, VX, VY, L, I, ds, alpha)

% POSTURE2RBM converts skeleton information and a drag coefficient ratio
% into a matrix of rigid body motions.

% Inputs
%   TX, TY     - x and y componenets of the unit tangent vectors. Matrix of
%                size (num skel points) x (number of frames)
%   DX,DY      - relative position of skeleton points with respect to the
%                center of mass
%   VX, VY     - velocities on the skeleton points without rigid body 
%                motion
%   L          - the total length of the skeleton
%   I          - moment of inertia of the skeleton. Vector of length number
%                of frames
%   ds         - the arclength increment between skeleton points
%   alpha      - the friction coefficient ratio (alpha_normal/alpha_tangent)
% 
% Outputs
%   RBM        - a matrix of rigid body motion. RBM(:,1:2) are x and y
%                components of translational velocity. RBM(:, 3) is the
%                angular velocity.  (number of frames - 1) x 3 matrix.
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
RBM = zeros(size(TX, 1)-1,3);

% get tangential component of velocity at each skeleton point without
% centre of mass (the wiggles)
TdotV = TX(2:end,:).*VX + TY(2:end,:).*VY;

% get cross product of relative skeleton position with tangent
DelXcrossT = DX(2:end,:).*TY(2:end,:) - DY(2:end,:).*TX(2:end,:);

for ii = 1:size(TX, 1)-1
    % assemble right hand side
    b1 = (alpha - 1).*ds.*trapz(TX(ii+1,:).*TdotV(ii,:));
    b2 = (alpha - 1).*ds.*trapz(TY(ii+1,:).*TdotV(ii,:));
    b3 = (alpha - 1).*ds.*trapz(DelXcrossT(ii,:).*TdotV(ii,:));
    
    % the matrix relating rigid body motion to the wiggles
    A11 = alpha*L + (1 - alpha).*ds.*trapz(TX(ii+1,:).^2);
    A12 = (1 - alpha).*ds.*trapz(TX(ii+1,:).*TY(ii+1,:));
    A13 = (1 - alpha).*ds.*trapz(TX(ii+1,:).*DelXcrossT(ii,:));
    
    A22 = alpha*L + (1 - alpha).*ds.*trapz(TY(ii+1,:).^2);
    A21 = A12;
    A23 = (1 - alpha).*ds.*trapz(TY(ii+1,:).*DelXcrossT(ii,:));
    
    A31 = A13;
    A32 = A23;
    A33 = alpha*I(ii+1) + (1 - alpha).*ds.*trapz(DelXcrossT(ii,:).^2);
    
    % solve the linear system
    bvec = [b1 b2 b3]';
    Amat = [A11 A12 A13; A21 A22 A23; A31 A32 A33];
    RBM(ii,:) = Amat'\bvec;
end