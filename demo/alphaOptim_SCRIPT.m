% This program computes the rigid body motion (translation and rotation) of
% a crawling worm from posture data.  The RBM is computed using a zero
% force and torque condition and linear resistive force model.

% ------------------------Set parameters-----------------------------------
L = 1.0; %worm length
dt = 1.0; %time between frames
alpha0 = 50; % initial guess of alpha
minFrameNum = 1000; % movies with fewer than this many frames are skipped
frameNum = 3000; % analyse only 1:frameNum frames in each movie
% -------------------------------------------------------------------------

% seed for reproducibility
rng(341)

% set the root directory
% directory = '/Users/abrown/Andre/wormVideos/results-12-05-10/Laura Grundy/';
directory = '/Users/abrown/Andre/wormVideos/results-12-05-10/wild-isolates/';

% get the file names
[fileList, ~] = dirSearch(directory, '_features.mat');

% intialise
alphaVec = NaN(numel(fileList), 1);
SDmat = NaN(numel(fileList), frameNum-1);
speedMat = NaN(numel(fileList), frameNum);
bendMat = NaN(numel(fileList), frameNum);


% loop through files
for jj = 1:numel(fileList)
    disp(jj/numel(fileList))

    % ------------------------ Load/prepare data --------------------------

    load(fileList{jj})
    
    % skip file if it has too few points
    if size(worm.posture.skeleton.x, 2) < frameNum
        continue
    end
    if sum(~isnan(worm.posture.skeleton.x(1, 1:frameNum))) < minFrameNum
        continue
    end
    
    % save speed and bend data
    speedMat(jj, :) = worm.locomotion.velocity.midbody.speed(1:frameNum);
    bendMat(jj, :) = worm.posture.bends.midbody.mean(1:frameNum);

    % calculate mean worm length
    mean_worm_L = nanmean(worm.morphology.length);
    
    % pre-process skeletons
    [X, Y] = preprocSkel(worm.posture.skeleton.x, ...
        worm.posture.skeleton.y, mean_worm_L);
    
    % check if preprocessing reduced number of frames too much
    if size(X, 1) < frameNum
        speedMat(jj, :) = NaN(size(speedMat(jj, :)));
        bendMat(jj, :) = NaN(size(bendMat(jj, :)));
        continue
    end
    
    X = X(1:frameNum, :);
    Y = Y(1:frameNum, :);

    % calculate arclength increment
    ds = L/(size(X, 2)-1);

    
    % calculate the observed rigid body motion
    [XCM, YCM, UX, UY, UXCM, UYCM, TX, TY, NX, NY, I, OMEG] = ...
        getRBM(X, Y, L, ds, dt);
    
    % subtract the observed rigid body motion
    [DX, DY, ODX, ODY, VX, VY, Xtil, Ytil, THETA] = ...
        subtractRBM(X, Y, XCM, YCM, UX, UY, UXCM, UYCM, OMEG, dt);
    
    % ---------------------- Running the model ----------------------------
    
    alphaVec(jj) = fminsearch(@(alpha) ...
        RBM_SDseg(UX, UY, TX, TY, DX, DY, VX, VY, L, I, ds, alpha), ...
        alpha0, optimset('TolFun', 1e-2, 'TolX', 1e-2));
    
    % use the optimal alpha to predict the rigid body motion
    RBM = posture2RBM(TX, TY, DX, DY, VX, VY, 1, I, ds, alphaVec(jj));
    
    % calculate error
    SD = SDseg(UX, UY, DX, DY, VX, VY, RBM);
    SDmat(jj, :) = SD';

end % loop over files
