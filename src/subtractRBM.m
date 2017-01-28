function [DX, DY, ODX, ODY, VX, VY, Xtil, Ytil, THETA] = ...
    subtractRBM(X, Y, XCM, YCM, UX, UY, UXCM, UYCM, OMEG, dt)

% SUBTRACTRBM takes two matrices of x and y coordinates and returns an
% equivalent set of coordinates with the rigid body motion (translation and
% rotation) subtracted away, leaving just the posture information.

% Inputs
%   X,Y - the x and y coordinates of the skeleton
%   XCM, YCM - the x and y coordinates of the centre of mass of the
%              skeleton.  Vectors with length of frame number.
%   UX, UY   - velocity components of the skeleton points.  Matrix of
%              (num skel points) x (number of frames - 1)
%   UXCM, UYCM - the x and y components of the centre of mass velocity
%   OMEG       - angular velocity of the skeleton. Vector of (length number
%                of frames - 1)
%   dt         - the time between frames
% 
% Outputs
%   DX,DY      - relative position of skeleton points with respect to the
%                center of mass
%   ODX, ODY   - velocity of skeleton points due to angular velocity. Note
%                sign of ODX is negative (taken into account in calculation
%                of VY)
%   VX, VY     - velocities on the skeleton points without rigid body 
%                motion
%   Xtil, Ytil - Rotated values of DX and DY.
%   THETA       - Rotation angle used to obtain Xtil and Ytil
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
DX = zeros(size(X));
DY = zeros(size(X));
ODX = zeros(size(X, 1)-1,size(X, 2));
ODY = zeros(size(X, 1)-1,size(X, 2));

Xtil = zeros(size(X));
Ytil = zeros(size(X));
THETA = zeros(size(X,1),1);
THETA(1) = 0;

%Subtracting the rigid body motion (RBM) from the data
for ii = 1:size(X, 1)
    DX(ii,:) = X(ii,:) - XCM(ii);
    DY(ii,:) = Y(ii,:) - YCM(ii);
    Xtil(ii,:) = DX(ii,:);
    Ytil(ii,:) = DY(ii,:);
        
    if(ii > 1)
        % cross product of dX with U (for angular velocity)
        itm1 = ii-1;
        
        ODX(itm1,:) = OMEG(itm1).*DX(ii,:);
        ODY(itm1,:) = OMEG(itm1).*DY(ii,:);
        THETA(ii) = THETA(ii-1) + OMEG(itm1)*dt;
   
        Xtil(ii,:) = cos(THETA(ii))*DX(ii,:) + sin(THETA(ii))*DY(ii,:);
        Ytil(ii,:) = cos(THETA(ii))*DY(ii,:) - sin(THETA(ii))*DX(ii,:);
       
    end
   
end


VX = UX - repmat(UXCM,1,size(X,2)) + ODY;
VY = UY - repmat(UYCM,1,size(X,2)) - ODX;
