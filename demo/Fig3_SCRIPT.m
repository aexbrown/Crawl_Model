% script to compare the coefficient ratio and fit quality of the linear
% model across multiple strains.  Reproduces Fig. 3.


% load the wild isolate workspace
wild = load('/Users/abrown/Andre/wormVideos/eric/RabetsEtAlModel/alphaOptim_wildisolates_WORKSPACE.mat');

% exclude AQ2947 data (which is just N2)
dropInds = [];
for ii = 1:numel(wild.fileList)
    if ~isempty(strfind(wild.fileList{ii}, 'AQ2947'))
        dropInds = [dropInds, ii];
    end
end
wild.fileList(dropInds) = [];
wild.alphaVec(dropInds) = [];

% get the worm names
wormNamesWild = cell(numel(wild.fileList), 1);
for ii = 1:numel(wild.fileList)
    % get the positions of forward slashes in file name
    slashPositions = strfind(wild.fileList{ii}, '/');
    
    % get the strain name
    wormNamesWild{ii} = wild.fileList{ii}(slashPositions(7)+1:slashPositions(8)-1);
end


% load the mutant workspace
mut = load('/Users/abrown/Andre/wormVideos/eric/RabetsEtAlModel/alphaOptim_mutants_WORKSPACE.mat');

% exclude "no_wait" data and nRHO data
dropInds = [];
for ii = 1:numel(mut.fileList)
    if ~isempty(strfind(mut.fileList{ii}, 'no_wait'))
        dropInds = [dropInds, ii];
    end
    if ~isempty(strfind(mut.fileList{ii}, 'RHO'))
        dropInds = [dropInds, ii];
    end
    if ~isempty(strfind(mut.fileList{ii}, 'nca'))
        dropInds = [dropInds, ii];
    end
    if ~isempty(strfind(mut.fileList{ii}, ';'))
        dropInds = [dropInds, ii];
    end
end
mut.fileList(dropInds) = [];
mut.alphaVec(dropInds) = [];

% get the worm names
wormNamesMut = cell(numel(mut.fileList), 1);
for ii = 1:numel(mut.fileList)
    % check for N2 hermaphrodites
    if ~isempty(strfind(mut.fileList{ii}, '/N2/'))
        if ~isempty(strfind(mut.fileList{ii}, '/XO/'))
            wormNamesMut{ii} = 'N2_male';
        else
            wormNamesMut{ii} = 'N2_herm';
        end
        continue
    end
        
    % get the positions of forward slashes in file name
    slashPositions = strfind(mut.fileList{ii}, '/');
    
    % get the string with the strain or mutant/allele name
    wormName = mut.fileList{ii}(slashPositions(7)+1:slashPositions(9)-1);
    
    %  remove '/on_food/' (present in wild isolate names and replace 
    % slashes with underscores
    wormName = strrep(wormName, '/on_food', '');
    wormName = strrep(wormName, '/', '_');
    wormNamesMut{ii} = wormName;
end


% join the wild isolate and mutant worm names and alpha vectors
wormNames = [wormNamesWild; wormNamesMut];
alphaVec = [wild.alphaVec; mut.alphaVec];
speedVec = [nanmean(abs(wild.speedMat), 2); nanmean(abs(mut.speedMat), 2)];
bendVec = [nanmean(abs(wild.bendMat), 2); nanmean(abs(mut.bendMat), 2)];
rmsdVec = [sqrt(nanmean(wild.SDmat, 2)); sqrt(nanmean(mut.SDmat, 2))];

% make vector to colour bars by either wild isolates or uncs
colorInds = zeros(numel(wormNames), 1);
colorInds(1:numel(wormNamesWild)) = 1;
uncMutantList = {'unc', 'acr-2_', 'egl-10_', 'egl-18_', 'egl-21_', ...
    'egl-2_', 'egl-31_', 'egl-33_', 'egl-47_', 'egl-5_', 'egl-5_', ...
    'egl-6_', 'gpb-2_'};
for ii = 1:numel(wormNames)
    for jj = 1:numel(uncMutantList)
        if ~isempty(strfind(wormNames{ii}, uncMutantList{jj}))
            colorInds(ii) = 2;
        end
    end
end
uniqueNames = unique(wormNames);

strainMedian = NaN(numel(uniqueNames), 1);
for ii = 1:numel(uniqueNames)
    % get the current inds
    currentInds = strcmp(wormNames, uniqueNames{ii});
    currentAlphas = alphaVec(currentInds);
    strainMedian(ii) = nanmedian(currentAlphas(currentAlphas < 200));
end
[sortedMedians, sortInds] = sort(strainMedian, 'descend');
namesOrdered = uniqueNames(sortInds);

% make a box plot of coefficient ratios
figure
boxplot(alphaVec(alphaVec < 200), wormNames(alphaVec < 200), ...
    'GroupOrder', namesOrdered, 'PlotStyle', 'compact', ...
    'MedianStyle', 'line', 'ColorGroup', colorInds(alphaVec < 200));
h = findobj(gca,'tag','Outliers');
set(h, 'Visible', 'off')
ylim([0, 200])


% make box plots of speed, bend, and root mean square deviation using the 
% same ordering and colouring
figure
boxplot(speedVec(alphaVec < 200), wormNames(alphaVec < 200), ...
    'GroupOrder', namesOrdered, 'PlotStyle', 'compact', ...
    'MedianStyle', 'line', 'ColorGroup', colorInds(alphaVec < 200));
ylim([-20, 350])
h = findobj(gca,'tag','Outliers');
set(h, 'Visible', 'off')

figure
boxplot(bendVec(alphaVec < 200), wormNames(alphaVec < 200), ...
    'GroupOrder', namesOrdered, 'PlotStyle', 'compact', ...
    'MedianStyle', 'line', 'ColorGroup', colorInds(alphaVec < 200));
ylim([0, 40])
h = findobj(gca,'tag','Outliers');
set(h, 'Visible', 'off')

figure
boxplot(rmsdVec(alphaVec < 200), wormNames(alphaVec < 200), ...
    'GroupOrder', namesOrdered, 'PlotStyle', 'compact', ...
    'MedianStyle', 'line', 'ColorGroup', colorInds(alphaVec < 200));
ylim([0, 0.03])
h = findobj(gca,'tag','Outliers');
set(h, 'Visible', 'off')


