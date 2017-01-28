% script to plot several short trajectories and compare the linear model
% predictions for several values of alpha (ratio of normal to tangential
% friction coefficients).  Choose example trajectories that show forward
% locomotion, reversals, and sharp turns.

% set parameters
trajNum = 25; % how many trajectories of each type should be plotted
trajLength = 1000; % trajectory length in frames
forwardDistThresh = 1e3; % distance threshold to define 'forward' segments
speedThresh = 250; % speed threshold for motion segments
backDistThresh = 1e3; % distance threshold to define 'reversal' segments
alphaLimits = [7, 1e6]; % alpha values to use (in addition to bestfit value)


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
        if ~isempty(worm.locomotion.motion.forward.frames)
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
    end
    
    % if still need more, check for reversal segments
    if backNum < trajNum
        % check for backward motion segments
        if ~isempty(worm.locomotion.motion.backward.frames)
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
                    if nanmean(abs(speed(round(startInd - trajLength/2):round(startInd + trajLength/2)))) > ...
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
alphaOptimVec = NaN(size(indMat, 1), 1);
distVecNoSlip = NaN(size(indMat, 1), 1);
distVecExp = NaN(size(indMat, 1), 1);
distVecBestFit = NaN(size(indMat, 1), 1);
for ii = 1:size(indMat, 1)
    disp(ii/size(indMat, 1))
    % load the feature file
    load(fileList{indMat(ii, 1)})
    
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
    
    % use model to predict rigid body motion
    RBM = posture2RBM(TX, TY, DX, DY, VX, VY, 1, I, ds, alphaLimits(1));
    
    % add the predicted rigid body motion back to the shifted skeletons
    [XCMrecon, YCMrecon] = RBM2traj(RBM, 1);

    % get the distance from the small alpha case
    distVecExp(ii) = sum(sqrt(diff(XCMrecon).^2 + diff(YCMrecon).^2));
    
    % plot the reconstructed trajectory
    plot(XCMrecon, YCMrecon)
    
    % use model to predict rigid body motion
    RBM = posture2RBM(TX, TY, DX, DY, VX, VY, 1, I, ds, alphaLimits(2));
    
    % add the predicted rigid body motion back to the shifted skeletons
    [XCMrecon, YCMrecon] = RBM2traj(RBM, 1);
    
    % get the distance from the large alpha case
    distVecNoSlip(ii) = sum(sqrt(diff(XCMrecon).^2 + diff(YCMrecon).^2));
    
    % plot the reconstructed trajectory
    plot(XCMrecon, YCMrecon)
    
    % use model to predict rigid body motion
    % this time find optimal alpha using fminsearch
    alpha0 = 50; % initial guess of alpha
    alphaOptim = fminsearch(@(alpha) ...
        RBM_SDseg(UX, UY, TX, TY, DX, DY, VX, VY, 1, I, ds, alpha), ...
        alpha0, optimset('TolFun', 1e-2, 'TolX', 1e-2));
    alphaOptimVec(ii) = alphaOptim;
    
    RBM = posture2RBM(TX, TY, DX, DY, VX, VY, 1, I, ds, alphaOptim);
    
    % add the predicted rigid body motion back to the shifted skeletons
    [XCMrecon, YCMrecon] = RBM2traj(RBM, 1);

    % get the distance from the best-fit alpha case
    distVecBestFit(ii) = sum(sqrt(diff(XCMrecon).^2 + diff(YCMrecon).^2));
    
    % plot the reconstructed trajectory
    plot(XCMrecon, YCMrecon)
    
    % add a point to indicate location of event if on reversal or turn data
    if ii > trajNum
        plot(XCM(round(trajLength/2)) - XCM(1), YCM(round(trajLength/2)) - YCM(1), ...
            '.', 'Color', 'r', 'MarkerSize', 15)
    end
    
    axis off
    saveas(gcf, ['trajComparePlot_largeN_' num2str(ii) '.eps'], 'epsc')
    close
end


% make box plot of the relative distances
figure
relDistData = [distVecExp./distVecBestFit, distVecNoSlip./distVecBestFit];
boxplot(relDistData, {'Experimental Ratio'; 'No-slip'})


% make box plot of the optimal friction coefficient ratios in each
% condition
figure
coeffData = [alphaOptimVec(1:trajNum), ...
    alphaOptimVec(trajNum+1:2*trajNum), alphaOptimVec(2*trajNum+1:end)];
coeffData(coeffData > 200) = NaN;
boxplot(coeffData, {'Forward'; 'Reversals'; 'Turns'});
p1 = ranksum(coeffData(:, 1), coeffData(:, 2));
p2 = ranksum(coeffData(:, 1), coeffData(:, 3));
p3 = ranksum(coeffData(:, 2), coeffData(:, 3));
disp([p1, p2, p3]*3) % Bonferroni adjusted p-vals. For n = 25: 0.1426    0.1968    2.8325
