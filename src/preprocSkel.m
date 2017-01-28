function [Xproc, Yproc] = preprocSkel(X, Y, meanL)

% PREPROCSKEL performs preprocessing on skeleton X and Y coordinate
% matrices.  Points are resampled to have more uniform arclength spacing,
% normalised by the mean length, and missing (NaN) data is filled in usng
% linear interpolation.  Matrices are also transposed.
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


% drop leading and trailing NaN values (don't want to extrapolate to
% extreme values)
if isnan(X(1, 1))
    % get the end of the starting NaN segment
    nanEnd = find(~isnan(X(1, :)), 1, 'first') - 1;
    
    % drop these values
    X(:, 1:nanEnd) = [];
    Y(:, 1:nanEnd) = [];
end

% is the last point NaN?
if isnan(X(1, end))
    % get the start of the final NaN segment
    nanStart = find(~isnan(X(1, :)), 1, 'last') + 1;
    
    % drop these values
    X(:, nanStart:end) = [];
    Y(:, nanStart:end) = [];
end

% reinterpolate X and Y to have more uniform spacing in terms of
% arclength
Xi = NaN(size(X));
Yi = NaN(size(Y));
for ii = 1:size(X, 2)
    if ~isnan(X(1, ii))
        t = cumsum(sqrt([0,diff(X(:, ii)')].^2 + [0,diff(Y(:, ii)')].^2));
        ti = linspace(t(1),t(end),size(X, 1));
        Xi(:, ii) = interp1(t,X(:, ii)',ti,'PCHIP');
        Yi(:, ii) = interp1(t,Y(:, ii)',ti,'PCHIP');
    end
end

Xi = Xi./meanL;
Yi = Yi./meanL;

%centreline interpolation for missing data
Xproc = NoNaNs(Xi)';
Yproc = NoNaNs(Yi)';