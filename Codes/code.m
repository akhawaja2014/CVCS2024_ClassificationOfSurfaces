

% PATHS for real data
testimage_path = 'E:/gracqrti/commodus-stela-w946/Experiments/MLIC/0001.png'
%RTI_acq = 'E:/gracqrti/Ex5/it6/acq'
RTI_acq = 'E:/gracqrti/commodus-stela-w946/Experiments/MLIC'
%figurefolder = 'E:/gracqrti/data/Experiments/Ex7/it2/figs'
%figurepath = 'E:/gracqrti/data/Experiments/Ex7'
lpfile_path = 'E:/gracqrti/Manuscript/lp27.lp'




%%

I = imread(testimage_path);  % Read the image file
imshow(I);  % Display the image
h = imrect;  % Interactively draw a rectangle
% Get the position of the rectangle
roiPosition = wait(h);
% Create a binary mask for the ROI
mask = false(size(I, 1), size(I, 2));
mask(round(roiPosition(2)):round(roiPosition(2)+roiPosition(4)), round(roiPosition(1)):round(roiPosition(1)+roiPosition(3))) = true;
% Display the mask (optional)
figure;
imshow(mask);

%imae = imread(testimage_path);
% Display the test image
%figure;
%imshow(imae);
%impixelinfo;

%% -----------------------------------
% Read the RTI data folder
folder = RTI_acq; % Change the path of folder. 
filePattern = fullfile(folder, '*.png');  % Just to detect images in the folder.
theFiles = dir(filePattern); % lists files and folders in the current folder

mlic = {}; % Initialize cell array for reading all the RTI images 

for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name; %theFiles is a structure array whose field 'name'(e.g. 0001.png) is accessed here
    fullFileName = fullfile(folder, baseFileName); % append the path
    fprintf(1, 'Now reading %s\n', fullFileName); % just to tell that its reading image, everything okay till now.
    img_clr = imread(fullFileName); % read image
    %img = im2gray(img_clr); % convert to grayscale if not
    if ~isempty(img_clr)
        mlic{end+1} = img_clr; % if the variable img is not empty append the cell array mlic 
    end
end

% How to detect most dynamic pixels?

%% ----------------------------------
%get pixel location of mask
% Assuming 'mask' is your binary mask...
tic
[row, col] = find(mask);

% Initialize a cell array to store ROI signals
ROI_signals = cell(1, length(row));
for idx = 1:length(row)
    % Get the current pixel location
    currRow = row(idx);
    currCol = col(idx);
    
    % Initialize an array to store the pixel intensities
    pixelIntensities = zeros(1, length(mlic));
    
    for i = 1 : length(mlic)
        pixelIntensities(i) = mlic{i}(currRow, currCol);  % MATLAB indexing is 1-based not 0-based
    end
    
    % Store the pixel intensities in the ROI_signals cell array
    ROI_signals{idx} = pixelIntensities;
end
% Read the signal in RTI stack of images

% Calculate the standard deviation of each signal in ROI_signals
stdDevs = cellfun(@std, ROI_signals);



% Assume ROI_signals is your data and each cell contains one signal

%%

numSignals = length(ROI_signals);

% Normalize and ensure all signals have the same length
maxLength = max(cellfun(@length, ROI_signals));
dataMatrix = zeros(numSignals, maxLength);

for i = 1:numSignals
    currentSignal = ROI_signals{i};
    dataMatrix(i, 1:length(currentSignal)) = interp1(1:length(currentSignal), currentSignal, linspace(1, length(currentSignal), maxLength), 'linear');
end

% Normalize each signal to range [0,1]
%dataMatrix = (dataMatrix - min(dataMatrix, [], 2)) ./ (max(dataMatrix, [], 2) - min(dataMatrix, [], 2));

% Define SOM parameters
mapDimensions = [2, 2]; % e.g., 10x10 grid. Adjust based on your needs
epochs = 100; % Number of training iterations

% Train SOM
som_net = selforgmap(mapDimensions, epochs);
som_net = train(som_net, dataMatrix');

% Get clusters
outputs = som_net(dataMatrix');
clusters = vec2ind(outputs);

%% Plotting the SOM signals toghether in o

% % Colors representing each cluster
% colors = {'r', 'g', 'b', 'y', 'm', 'c', 'k'}; % extend this if you have more clusters
% 
% figure;
% 
% % Get the unique cluster assignments
% uniqueClusters = unique(clusters);
% 
% % To store handles for legend entries
% legendHandles = zeros(length(uniqueClusters), 1);
% 
% % Loop through each unique cluster
% for i = 1:length(uniqueClusters)
%     % Extract the signals that belong to the current cluster
%     currentClusterSignals = ROI_signals(clusters == uniqueClusters(i));
%     
%     % Plot each signal in the current cluster
%     for j = 1:length(currentClusterSignals)
%         plot(currentClusterSignals{j}, 'Color', colors{mod(i-1, length(colors)) + 1}, 'LineWidth', 0.5);
%         hold on;
%     end
%     
%     % Plot a dummy line for legend (invisible line with a label for the legend)
%     legendHandles(i) = plot(NaN,NaN, 'Color', colors{mod(i-1, length(colors)) + 1}, 'LineWidth', 2);
% end
% 
% legend(legendHandles, arrayfun(@(x) ['Cluster ', num2str(x)], uniqueClusters, 'UniformOutput', false));
% hold off;
% 
% 
% % Assuming 'clusters' contains the class indices for each signal
% number_Signals = length(clusters);


%% Plotting SOM clusters seperately

% Colors representing each cluster
colors = {'r', 'g', 'b', 'y', 'm', 'c', 'k'}; % extend this if you have more clusters

% Get the unique cluster assignments
uniqueClusters = unique(clusters);

% Loop through each unique cluster
for i = 1:length(uniqueClusters)
    figure; % Create a new figure for each cluster
    
    % Extract the signals that belong to the current cluster
    currentClusterSignals = ROI_signals(clusters == uniqueClusters(i));
    
    % Plot each signal in the current cluster
    for j = 1:length(currentClusterSignals)
        plot(currentClusterSignals{j}, 'Color', colors{mod(i-1, length(colors))+1}, 'LineWidth', 0.1);
        hold on;
    end
    
    % Plot a dummy line for legend (invisible line with a label for the legend)
    legendHandle = plot(NaN, NaN, 'Color', colors{mod(i-1, length(colors)) + 1}, 'LineWidth', 2);
    
    % Add the legend
    legend(legendHandle, ['Cluster ', num2str(uniqueClusters(i))]);
    % Improve the appearance of the axes
    ax = gca;
    ax.FontSize = 15;
    ax.LineWidth = 2;
    ax.Box = 'on'; % Add a box around the plot
    % Set the x-axis limits
    xlabel('MLIC Image no');
    ylabel('Signal Intensity');
    xlim([1 27]);
    ylim([0 255]);
    % Add grid for better readability
    grid on;
    hold off;
end



% Assuming 'clusters' contains the class indices for each signal
number_Signals = length(clusters);




%%


% Use k-means clustering to cluster the signals based on their standard deviations
numClusters = 4;
stdDevs = stdDevs(:); % Ensure stdDevs is a column vector


figure; % Create a new figure window

% Plot the standard deviations using scatter plot
scatter(1:length(stdDevs), stdDevs);

title('Standard Deviation of Signals from ROI'); % Set the title
xlabel('Signal Index'); % X-axis label
ylabel('Standard Deviation'); % Y-axis label
grid on; % Display grid lines for better visualization


options = statset('MaxIter', 5000000, 'TolFun', 1e-7);
[cluster_idx, ~] = kmeans(stdDevs, numClusters, 'Options', options);


% Now, the 'cluster_idx' contains the cluster assignment for each signal.
% cluster_idx = 1 indicates the first cluster and cluster_idx = 2 indicates the second cluster.

% If you want to see which signals belong to the first cluster:
signals_cluster1 = ROI_signals(cluster_idx == 1);

% If you want to see which signals belong to the second cluster:
signals_cluster2 = ROI_signals(cluster_idx == 2);

%%
% % Find the indices of signals that belong to cluster 1
% cluster1_indices = find(cluster_idx == 1);
% 
% % Extract signals belonging to cluster 1
% cluster1_signals = ROI_signals(cluster1_indices);
% 
% % Plot signals from cluster 1
% figure;
% hold on;
% for i = 1:length(cluster1_signals)
%     plot(cluster1_signals{i});
% end
% hold off;
% title('Signals belonging to Cluster 1');
% xlabel('MLIC Image no');
% ylabel('Signal Intensity');
% 
% 
% %%
% 
% % Find the indices of signals that belong to cluster 1
% cluster2_indices = find(cluster_idx == 2);
% 
% % Extract signals belonging to cluster 1
% cluster2_signals = ROI_signals(cluster2_indices);
% 
% % Plot signals from cluster 1
% figure;
% hold on;
% for i = 1:length(cluster2_signals)
%     plot(cluster2_signals{i});
% end
% hold off;
% title('Signals belonging to Cluster 2');
% xlabel('MLIC Image no');
% ylabel('Signal Intensity');
% 
% %%
% % Find the indices of signals that belong to cluster 1
% cluster3_indices = find(cluster_idx == 3);
% 
% % Extract signals belonging to cluster 1
% cluster3_signals = ROI_signals(cluster3_indices);
% 
% % Plot signals from cluster 1
% figure;
% hold on;
% for i = 1:length(cluster3_signals)
%     plot(cluster3_signals{i});
% end
% hold off;
% title('Signals belonging to Cluster 3');
% xlabel('MLIC Image no');
% ylabel('Signal Intensity');
% 
% %%
% % Find the indices of signals that belong to cluster 1
% cluster4_indices = find(cluster_idx == 4);
% 
% % Extract signals belonging to cluster 1
% cluster4_signals = ROI_signals(cluster4_indices);
% 
% % Plot signals from cluster 1
% figure;
% hold on;
% for i = 1:length(cluster4_signals)
%     plot(cluster4_signals{i});
% end
% hold off;
% title('Signals belonging to Cluster 4');
% xlabel('MLIC Image no');
% ylabel('Signal Intensity');

%%  k means plot 

% Define the number of clusters
num_clusters = 4;

% Colors representing each cluster
colors = {'r', 'g', 'b', 'y', 'm', 'c', 'k'}; % extend this if you have more clusters

% Loop through each cluster and plot the signals
for cluster_id = 1:num_clusters
    % Find the indices of signals that belong to the current cluster
    cluster_indices = find(cluster_idx == cluster_id);

    % Extract signals belonging to the current cluster
    cluster_signals = ROI_signals(cluster_indices);

    % Plot signals from the current cluster
    figure;
    hold on;
    for i = 1:length(cluster_signals)
        plot(cluster_signals{i}, 'Color', colors{cluster_id});
    end
    hold off;
    
    % Set plot title and labels
    %title(['Signals belonging to Cluster ', num2str(cluster_id)]);
    legend(['Cluster ', num2str(cluster_id)]);
    xlabel('MLIC Image no');
    ylabel('Signal Intensity');
    ax = gca;
    ax.FontSize = 15;
    ax.LineWidth = 2;
    ax.Box = 'on'; % Add a box around the plot
    % Set the x-axis limits
    xlim([1 27]);
    ylim([0 255]);
    % Add grid for better readability
    grid on;

    
end




%%    SOM


% Assuming clusters, row, and col are already defined

% Create a transparent mask based on the clusters
mask_img = zeros(size(I, 1), size(I, 2), 3);  % Initialize mask

% Define colors for each cluster (R, G, B)
colors = [
    1, 0, 0;      % Red
    0, 1, 0;      % Green
    0, 0, 1;      % Blue
    1, 1, 0;      % Yellow
    0, 1, 1;      % Cyan
    1, 0, 1;      % Magenta
    0.5, 0, 0;    % Dark Red
    0, 0.5, 0;    % Dark Green
    0, 0, 0.5;    % Dark Blue
    0.5, 0.5, 0;  % Olive
    0, 0.5, 0.5;  % Teal
    0.5, 0, 0.5;  % Purple
    0.7, 0.7, 0.7; % Grey
    0.8, 0.5, 0;  % Orange
    0.5, 0.8, 0.8; % Pale Cyan
    0.8, 0.8, 0.5; % Pale Yellow
    0.6, 0.2, 0.2; % Brown
    0.9, 0.6, 0.7; % Pink
    0, 0.7, 0.7;  % Turquoise
    0.7, 0.7, 1;  % Light Blue
];


% Assign colors to the mask based on clusters
for idx = 1:length(row)
    cluster_id = clusters(idx);
    mask_img(row(idx), col(idx), :) = colors(cluster_id, :);
end

% Display the original image
figure;
imshow(I);
hold on;

% Overlay the mask with transparency
h = imshow(mask_img);
set(h, 'AlphaData', 0.3);

% Title and show
title('Original Image with SOM Classification Masks','FontSize', 16);

hold off;


%% K Means

% Assuming clusters, row, and col are already defined

% Create a transparent mask based on the clusters
mask_img = zeros(size(I, 1), size(I, 2), 3);  % Initialize mask

% Define colors for each cluster (R, G, B)
colors = [
    1, 0, 0;      % Red
    0, 1, 0;      % Green
    0, 0, 1;      % Blue
    1, 1, 0;      % Yellow
    0, 1, 1;      % Cyan
    1, 0, 1;      % Magenta
    0.5, 0, 0;    % Dark Red
    0, 0.5, 0;    % Dark Green
    0, 0, 0.5;    % Dark Blue
    0.5, 0.5, 0;  % Olive
    0, 0.5, 0.5;  % Teal
    0.5, 0, 0.5;  % Purple
    0.7, 0.7, 0.7; % Grey
    0.8, 0.5, 0;  % Orange
    0.5, 0.8, 0.8; % Pale Cyan
    0.8, 0.8, 0.5; % Pale Yellow
    0.6, 0.2, 0.2; % Brown
    0.9, 0.6, 0.7; % Pink
    0, 0.7, 0.7;  % Turquoise
    0.7, 0.7, 1;  % Light Blue
];


% Assign colors to the mask based on clusters
for idx = 1:length(row)
    cluster_id = cluster_idx(idx);
    mask_img(row(idx), col(idx), :) = colors(cluster_id, :);
end

% Display the original image
figure;
imshow(I);
hold on;

% Overlay the mask with transparency
h = imshow(mask_img);
set(h, 'AlphaData', 0.3);

% Title and show
title('Original Image with K means Classification Masks','FontSize', 16);
hold off;
