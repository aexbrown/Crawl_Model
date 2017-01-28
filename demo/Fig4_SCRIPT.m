% script to plot several short trajectories and compare the nonlinear model
% predictions for several values of the model coefficients.  Choose example
% trajectories that show forward locomotion, reversals, and sharp turns.
% Reproduces Figure 4

% first find the appropriate value of the parameters by linear
% interpolation. Our data were collected with 2% agar plates.

agarPct = [1.0, 1.7, 3.0, 4.0, 6.0];
alphaT = [0.35, 0.32, 0.26, 0.24, 0.20];
alphaN = [4.6, 6.8, 7.9, 7.3, 5.8];
gammaT = [0.55, 0.58, 0.68, 0.77, 0.81];
gammaN = [0.27, 0.31, 0.30, 0.34, 0.38];

alphaT_2pct = interp1(agarPct, alphaT, 2); disp(alphaT_2pct);
alphaN_2pct = interp1(agarPct, alphaN, 2); disp(alphaN_2pct);
gammaT_2pct = interp1(agarPct, gammaT, 2); disp(gammaT_2pct);
gammaN_2pct = interp1(agarPct, gammaN, 2); disp(gammaN_2pct);

% set parameters
trajNum = 5; % how many trajectories of each type should be plotted
trajLength = 1000; % trajectory length in frames
forwardDistThresh = 1e3; % distance threshold to define 'forward' segments
speedThresh = 250; % speed threshold for motion segments
backDistThresh = 1e3; % distance threshold to define 'reversal' segments


% seed for reproducibility
rng(55551)

% set the root directory
directory = '/Users/abrown/Andre/wormVideos/results-12-05-10/wild-isolates/';

% get the file names
[fileList, ~] = dirSearch(directory, '_features.mat');


% get a random subset of the data
randInds = randi(numel(fileList), numel(fileList), 1);

% initialise
indMat = NaN(trajNum*3, 3);

forNum = 0;
backNum = 0;
turnNum = 0;

% loop through files to find appropriate segments
for ii = 1:numel(fileList)
    disp([ii, forNum, backNum, turnNum])
    % load a feature file
    load(fileList{randInds(ii)})
    
    % get the absolute speed for this worm
    speed = worm.locomotion.velocity.midbody.speed;
    
    % if still need more, check for forward motion segments
    if forNum < trajNum
        % check for forward motion segments
        forwardDistances = ...
            vertcat(worm.locomotion.motion.forward.frames.distance);
        matches = find(forwardDistances > forwardDistThresh);
        
        % check the candidates for speed criteria
        for jj = 1:numel(matches)
            % get current segment (add one because first frame seems to be
            % 0)
            startInd = ...
                worm.locomotion.motion.forward.frames(matches(jj)).start + 1;
            
            if startInd + trajLength - 1 <= length(speed)
                if nanmean(speed(startInd+1:startInd+trajLength)) > ...
                        speedThresh
                    % segment is a match
                    forNum = forNum + 1;
                    indMat(forNum, 1) = randInds(ii);
                    indMat(forNum, 2) = startInd;
                    indMat(forNum, 3) = startInd + trajLength - 1;
                    break
                end
            end
        end
    end
    
    % if still need more, check for reversal segments
    if backNum < trajNum
        % check for backward motion segments
        backwardDistances = ...
            vertcat(worm.locomotion.motion.backward.frames.distance);
        matches = find(backwardDistances > backDistThresh);
        
        % check the candidates for speed criteria
        for jj = 1:numel(matches)
            % get current segment
            startInd = worm.locomotion.motion.backward.frames(matches(jj)).start;
            
            % only continue if startInd is greater than half trajLength
            % because want reversal start centered in trajectory
            if round(startInd - trajLength/2) >= 1 ...
                    && round(startInd + trajLength/2) <= length(speed)
                if nanmean(abs(speed(round(startInd - trajLength/2):round(startInd + trajLength/2)))) > ...
                        speedThresh
                    % segment is a match
                    backNum = backNum + 1;
                    indMat(backNum+trajNum, 1) = randInds(ii);
                    % shift start and end so reversal start is roughly in
                    % middle of segment
                    indMat(backNum+trajNum, 2) = round(startInd - trajLength/2);
                    indMat(backNum+trajNum, 3) = round(startInd + trajLength/2);
                    break
                end
            end
        end
    end
    
    % if still need more, check for turn segments
    if turnNum < trajNum
        % check for turn motion segments
        if ~isempty(worm.locomotion.turns.omegas.frames)
            
            matches = vertcat(worm.locomotion.turns.omegas.frames.start);
            
            % check the candidates for speed criteria
            for jj = 1:numel(matches)
                % get current segment
                startInd = matches(jj);
                
                % only continue if startInd is greater than half trajLength
                % because want turn start centered in trajectory
                if round(startInd - trajLength/2) >= 1 ...
                    && round(startInd + trajLength/2) <= length(speed)
                    if nanmean(abs(speed(startInd+1:startInd+trajLength))) > ...
                            speedThresh
                        % segment is a match
                        turnNum = turnNum + 1;
                        indMat(turnNum+trajNum*2, 1) = randInds(ii);
                        % shift start and end so reversal start is roughly in
                        % middle of segment
                        indMat(turnNum+trajNum*2, 2) = round(startInd - trajLength/2);
                        indMat(turnNum+trajNum*2, 3) = round(startInd + trajLength/2);
                        break
                    end
                end
            end
        end
    end
    
    % check if enough segments have been found
    if forNum == trajNum && backNum == trajNum && turnNum == trajNum
        disp([ii, forNum, backNum, turnNum])
        break
    end
end


% plot the segments in each category
phiVec = NaN(size(indMat, 1), 1);
for ii = 1:size(indMat, 1)
    disp(ii/size(indMat, 1))
    % load the feature file
    load(fileList{indMat(ii, 1)})
    
    % use frame rate and worm length to calculate the parameter conversion 
    % from our unit length and time interval so that the parameters match 
    % those from Rabets et al.
    convFact = nanmean(worm.morphology.length)*info.video.resolution.fps;
    
    
    % first pre-process skeleton data
    [X, Y] = preprocSkel(...
        worm.posture.skeleton.x(:, indMat(ii, 2):indMat(ii, 3)), ...
        worm.posture.skeleton.y(:, indMat(ii, 2):indMat(ii, 3)), ...
        nanmean(worm.morphology.length));
    
    % calculate arclength increment
    ds = 1/(size(X, 2)-1);
    
    % calculate the observed rigid body motion
    [XCM, YCM, UX, UY, UXCM, UYCM, TX, TY, NX, NY, I, OMEG] = ...
        getRBM(X, Y, 1, ds, 1);
    
    % plot original centre of mass trajectory
    figure
    plot(XCM - XCM(1), YCM - YCM(1))
    xlim([-10, 10])
    ylim([-10, 10])
    axis equal
    hold on
    
    % subtract the observed rigid body motion
    [DX, DY, ODX, ODY, VX, VY, Xtil, Ytil, THETA] = ...
        subtractRBM(X, Y, XCM, YCM, UX, UY, UXCM, UYCM, OMEG, 1);
    
    % Calculate the tangential and normal velocities
    VT = VX.*TX(1:end-1,:) + VY.*TY(1:end-1,:);
    VN = VX.*NX(1:end-1,:) + VY.*NY(1:end-1,:);
    
    % use model to predict rigid body motion. Convert alpha paramters from
    % Rabets et al. using conversion factor
    RBM = posture2RBM_nonlinear(X, Y, XCM, YCM, TX, TY, VT, VN, NX, NY, ...
        ds, alphaT_2pct * convFact^gammaT_2pct, ...
        alphaN_2pct * convFact^gammaN_2pct, gammaT_2pct, gammaN_2pct);
    
    % add the predicted rigid body motion back to the shifted skeletons
    [XCMrecon, YCMrecon] = RBM2traj(RBM, 1);

    % plot the reconstructed trajectory
    plot(XCMrecon, YCMrecon)
    
    
    % use model to predict rigid body motion
    % this time find optimal parameters using fminsearch
    phi0 = 1.15; % initial guess of alpha
    phiOptim = fminsearch(@(phi) ...
        RBM_SD_nonlinear_phi(UXCM, UYCM, X, Y, XCM, YCM, TX, TY, ...
        VT, VN, NX, NY, ds, ...
        alphaT_2pct, alphaN_2pct, gammaT_2pct, gammaN_2pct, phi), ...
        phi0, optimset('TolFun', 1e-2, 'TolX', 1e-2));
    phiVec(ii) = phiOptim;
    
    % adjust parameters by phi
    [atS, anS, gtS, gnS] = phiTransform(alphaT_2pct, alphaN_2pct, ...
        gammaT_2pct, gammaN_2pct, phiOptim);

    RBM = posture2RBM_nonlinear(X, Y, XCM, YCM, TX, TY, VT, VN, NX, NY, ...
        ds, atS, anS, gtS, gnS);
    
    % add the predicted rigid body motion back to the shifted skeletons
    [XCMrecon, YCMrecon] = RBM2traj(RBM, 1);
    
    % plot the reconstructed trajectory
    plot(XCMrecon, YCMrecon)
    
    % add a point to indicate location of event if on reversal or turn data
    if ii > trajNum
        plot(XCM(round(trajLength/2)) - XCM(1), YCM(round(trajLength/2)) - YCM(1), ...
            '.', 'Color', 'r', 'MarkerSize', 15)
    end
    
    axis off
    saveas(gcf, ['trajComparePlot_' num2str(ii) '.eps'], 'epsc')
    close
end