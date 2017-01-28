function [XCM, YCM, UX, UY, UXCM, UYCM, TX, TY, NX, NY, I, OMEG] = ...
    getRBM(X, Y, L, ds, dt)

% GETRBM takes in matrices of x and y skeleton coordinates and calculates
% their rigid body motion (translation and rotation). It outputs the point
% velocities, tangent angles, and center of mass velocities.
% 
% Inputs
%   X,Y - the x and y coordinates of the skeleton
%   L   - the total length of the skeleton
%   ds  - the arclength increment between skeleton points
%   dt  - the time between frames
% 
% Outputs
%   XCM, YCM   - the x and y coordinates of the centre of mass of the
%                skeleton.  Vectors with length of frame number.
%   UX, UY     - velocity components of the skeleton points.  Matrix of
%                (num skel points) x (number of frames - 1)
%   UXCM, UYCM - the x and y components of the centre of mass velocity
%   TX, TY     - x and y componenets of the unit tangent vectors. Matrix of
%                size (num skel points) x (number of frames)
%   NX, NY     - x and y componenets of the unit normal vectors. Matrix of
%                size (num skel points) x (number of frames)
%   I          - moment of inertia of the skeleton. Vector of length number
%                of frames
%   OMEG       - angular velocity of the skeleton. Vector of (length number
%                of frames - 1)
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



%Center of mass position
XCM = ds*trapz(X,2)./L;
YCM = ds*trapz(Y,2)./L;

%Velocity of the centerline from the data
UX = (X(2:end,:) - X(1:end-1,:))/dt;
UY = (Y(2:end,:) - Y(1:end-1,:))/dt;

%Velocity of the center of mass from the data
UXCM = ds*trapz(UX,2)./L;
UYCM = ds*trapz(UY,2)./L;

%unit vector tangent to the worm
[TX, TY] = TvecFromSkel(X',Y',L,ds);

TX = -TX';
TY = -TY';
NX = -TY; %unit vector normal to the worm
NY = TX;

%Computing the angular velocity and moment of inertia from the data
Iint = zeros(size(X, 1),size(X, 2));
Omegint = zeros(size(X, 1),size(X, 2));

for ii = 1:size(X, 1)
    % get integrand for moment of intertia
    Iint(ii,:) = (X(ii,:) - XCM(ii)).^2 + (Y(ii,:) - YCM(ii)).^2;
    
    if(ii > 1)
        % cross product of dX with U (for angular velocity)
        itm1 = ii - 1;
        Omegint(itm1,:) = (X(ii,:) - XCM(ii)).*(UY(itm1,:)-UYCM(itm1)) - ...
            (Y(ii,:) - YCM(ii)).*(UX(itm1,:)-UXCM(itm1));
    end
end

I = ds*trapz(Iint,2); %moment of inertia
OMEG = ds*trapz(Omegint,2)./I; %angular velocity
OMEG = OMEG(1:end-1);
