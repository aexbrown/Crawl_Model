% script to reproduce Fig. 1 (schematic of model and illustration of
% procudure of subtracting and then calculating rigid body motion)


% load some worm data
load(['/Users/abrown/Andre/wormVideos/results-12-05-10/' ...
    'wild-isolates/ED3054/on_food/XX/30m_wait/L/tracker_1/' ...
    '2010-11-26___11_14_42/' ...
    '507 ED3054 on food R_2010_11_26__11_14_42___1___3_features.mat'])


% ------------------------------Part A-------------------------------------

% plot two skeletons to show motion
frame1 = 1205;
frame2 = 1230;
figure

% plot skeleton 1
plot(worm.posture.skeleton.x(:, frame1), worm.posture.skeleton.y(:, frame1), ...
    'LineWidth', 4)
hold on
plot(worm.posture.skeleton.x(1, frame1), worm.posture.skeleton.y(1, frame1), ...
    '.', 'MarkerSize', 15)

% plot skeleton 2
plot(worm.posture.skeleton.x(:, frame2), worm.posture.skeleton.y(:, frame2), ...
    'LineWidth', 4)
plot(worm.posture.skeleton.x(1, frame2), worm.posture.skeleton.y(1, frame2), ...
    '.', 'MarkerSize', 15)

% show lines connecting the equivalent points on both skeletons
skelMatX = [worm.posture.skeleton.x(:, frame1)'; ...
    worm.posture.skeleton.x(:, frame2)'];
skelMatY = [worm.posture.skeleton.y(:, frame1)'; ...
    worm.posture.skeleton.y(:, frame2)'];
plot(skelMatX, skelMatY, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5)

axis equal
hold off


% ------------------------------Part B-------------------------------------

% plot series of skeletons, subtract away rigid body motion, plot model
% result

% first pre-process skeleton data
[X, Y] = preprocSkel(worm.posture.skeleton.x, ...
    worm.posture.skeleton.y, nanmean(worm.morphology.length));

% calculate arclength increment (assumes constant arclength increment and
% total length 1)
ds = 1/(size(X, 2)-1);

% calculate the observed rigid body motion
[XCM, YCM, UX, UY, UXCM, UYCM, TX, TY, NX, NY, I, OMEG] = ...
    getRBM(X, Y, 1, ds, 1);

% subtract the observed rigid body motion
[DX, DY, ODX, ODY, VX, VY, Xtil, Ytil, THETA] = ...
    subtractRBM(X, Y, XCM, YCM, UX, UY, UXCM, UYCM, OMEG, 1);

% use model to predict rigid body motion
alpha = 86.1;
RBM = posture2RBM(TX, TY, DX, DY, VX, VY, 1, I, ds, alpha);

% add the predicted rigid body motion back to the shifted skeletons
[~, ~, Xrecon, Yrecon] = addRBM(DX, DY, VX, VY, RBM, 1);


% plot the skeletons in the different analysis stages
figure
plot(X(1:200, :)' - X(1, 1), Y(1:200, :)' - Y(1, 1), 'Color', 'r')
hold on
plot(X(1:200, 1)' - X(1, 1), Y(1:200, 1)' - Y(1, 1), ...
    '.', 'Color', 'b', 'MarkerSize', 15)
xlim([-0.8, 0.4])
ylim([-0.8, 0.6]) 
axis equal

figure
plot(DX(1:200, :)' - DX(1, 1), DY(1:200, :)' - DY(1, 1), 'Color', 'r')
hold on
plot(DX(1:200, 1)' - DX(1, 1), DY(1:200, 1)' - DY(1, 1), ...
    '.', 'Color', 'b', 'MarkerSize', 15)
xlim([-0.8, 0.4])
ylim([-0.8, 0.6]) 
axis equal

figure
plot(Xrecon(1:200, :)' - Xrecon(1, 1), Yrecon(1:200, :)' - Yrecon(1, 1), 'Color', 'r')
hold on
plot(Xrecon(1:200, 1)' - Xrecon(1, 1), Yrecon(1:200, 1)' - Yrecon(1, 1), ...
    '.', 'Color', 'b', 'MarkerSize', 15)
xlim([-0.8, 0.4])
ylim([-0.8, 0.6]) 
axis equal

