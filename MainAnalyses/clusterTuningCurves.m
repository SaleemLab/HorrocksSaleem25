%% ==================== PREPROCESSING & CLUSTERING ========================
% (Assumes you already have your data structure 'goodUnits' loaded.)

% --- Thresholding and Data Preparation ---
r2_thresh  = 0.1;
r2p_thresh = 0.05;

allTuning = cat(1, goodUnits.tuning_stat, goodUnits.tuning_run);
allr2     = [[goodUnits.r2_stat], [goodUnits.r2_run]]';
allr2p    = [[goodUnits.r2pval_stat], [goodUnits.r2pval_run]]';

% Select units that pass criteria:
idx         = find(allr2 >= r2_thresh & allr2p < r2p_thresh);
allCurves   = allTuning(idx, :);
allCurves_z = zscore(allCurves, [], 2);   % z-score each curve (row-wise)
nCurves     = size(allCurves_z, 1);

%% Sorting using optimalleaforder
% Compute Dissimilarity Matrix & Hierarchical Linkage
D  = pdist(allCurves_z, 'euclidean');
dm = squareform(D);
Z  = linkage(D, 'average');

% find the optimal leaf order
leafOrder = optimalleaforder(Z, D);

% re-order the dissimilarity matrix and tuning curves
dm_reordered = dm(leafOrder, leafOrder);
tuning_reordered = allCurves_z(leafOrder, :);



%% Clustering 

% use silhouette score to find numCluster
tic
clusterRange = 2:2:32
avgSilhouette = zeros(length(clusterRange),1);

for i = 1:length(clusterRange)
    % Compute cluster assignments for current number of clusters
    currentClusters = cluster(Z, 'maxclust', clusterRange(i));
    
    % Compute silhouette values for each observation
    s = silhouette(allCurves_z, currentClusters, 'euclidean');
    
    % Compute the average silhouette value
    avgSilhouette(i) = mean(s);
end

% Plot average silhouette value vs. number of clusters
figure;
plot(clusterRange, avgSilhouette, '-o', 'LineWidth', 2);
xlabel('Number of Clusters');
ylabel('Average Silhouette Value');
title('Silhouette Analysis for Optimal Number of Clusters');
grid on;
toc

%% --- Define number of clusters (change this variable as desired) ---
numClusters   = 8;   % <--- CHANGE THIS VALUE TO USE A DIFFERENT NUMBER OF CLUSTERS

% --- Obtain Cluster Assignments ---
cluster_idx = cluster(Z, 'maxclust', numClusters);
% Reorder the cluster assignments to match the leaf order:
clusters_reordered = cluster_idx(leafOrder);

% --- Determine boundaries in the sorted order ---
% These boundaries show where the cluster label changes
clusterBoundaries = find(diff(clusters_reordered)) + 0.5;

% --- Define a set of colors for the clusters (using a colormap with numClusters distinct colors) ---
clusterColors = lines(numClusters);  % numClusters-by-3 matrix of RGB values

% --- Determine the cutoff for the dendrogram ---
% One common choice is to take the merge distance at row (nCurves - numClusters)
cutoff = Z(nCurves - numClusters, 3);


%% ==================== FIGURE 1: DENDROGRAM & DISSIMILARITY MATRIX ====================
figure('Name', sprintf('Dendrogram & Reordered Dissimilarity Matrix (%d clusters)', numClusters), ...
       'Position', [100, 100, 1200, 600]);

% ----- Left Panel: Full Dendrogram with Cluster Coloring -----
subplot(1,2,1);
[H, T, outperm] = dendrogram(Z, nCurves, 'Reorder', leafOrder, 'Orientation', 'left', ...
    'ColorThreshold', cutoff);
title(sprintf('Dendrogram (%d final clusters)', numClusters));
xlabel('Distance');
set(gca, 'YTick', []);  % Remove y-axis tick labels for clarity
hold on;
% Draw a dashed vertical line at the cutoff
xline(cutoff, 'k--', 'LineWidth', 1);
hold off;

% ----- Right Panel: Reordered Dissimilarity Matrix with Cluster Boundaries -----
subplot(1,2,2);
dm_reordered = dm(leafOrder, leafOrder);
imagesc(dm_reordered);
axis square;
set(gca, 'YDir', 'normal');  % so that row 1 is at the top (matching the dendrogram)
colorbar;
title('Reordered Dissimilarity Matrix');
xlabel('Curve Index');
ylabel('Curve Index');
hold on;
% Draw white lines at the cluster boundaries (both horizontal and vertical)
for b = 1:length(clusterBoundaries)
    line([0.5, nCurves+0.5], [clusterBoundaries(b) clusterBoundaries(b)], 'Color', 'w', 'LineWidth', 2);
    line([clusterBoundaries(b) clusterBoundaries(b)], [0.5, nCurves+0.5], 'Color', 'w', 'LineWidth', 2);
end
hold off;

% ----- Add a Cluster Color Bar -----
% Create a narrow axes to display each curve's cluster label as a color bar.
axPos = get(gca, 'Position');  % position of current axes
cbWidth = 0.02;                % relative width for the color bar
cbPos = [axPos(1)+axPos(3)+0.01, axPos(2), cbWidth, axPos(4)];
axes('Position', cbPos);
imagesc(clusters_reordered(:)); % display as a column vector
colormap(gca, clusterColors);
set(gca, 'XTick', [], 'YTick', []);
title('Clusters');


%% ==================== FIGURE 2: DENDROGRAM & SORTED TUNING CURVES ====================
figure('Name', sprintf('Dendrogram & Sorted Tuning Curves (%d clusters)', numClusters), ...
       'Position', [150, 150, 1200, 600]);

% ----- Left Panel: Dendrogram (same as before) -----
subplot(1,2,1);
[H2, T2, outperm2] = dendrogram(Z, nCurves, 'Reorder', leafOrder, 'Orientation', 'left', ...
    'ColorThreshold', cutoff);
title(sprintf('Dendrogram (%d final clusters)', numClusters));
xlabel('Distance');
set(gca, 'YTick', []);
hold on;
xline(cutoff, 'k--', 'LineWidth', 1);
hold off;

% ----- Right Panel: Sorted (Z-scored) Tuning Curves with Cluster Boundaries -----
subplot(1,2,2);
tuning_reordered = allCurves_z(leafOrder, :);
imagesc(tuning_reordered);
axis tight;
set(gca, 'YDir', 'normal');  % so that row 1 is at the top, matching the dendrogram
colorbar;
title('Sorted Tuning Curves');
xlabel('Tuning Curve Bins');
ylabel('Curve Index');
hold on;
% Draw horizontal lines at the boundaries between clusters
for b = 1:length(clusterBoundaries)
    yline(clusterBoundaries(b), 'k-', 'LineWidth', 2);
end
hold off;

% ----- Add a Cluster Color Bar next to the Tuning Curves -----
axPos2 = get(gca, 'Position');
cbWidth = 0.02;
cbPos2 = [axPos2(1)+axPos2(3)+0.01, axPos2(2), cbWidth, axPos2(4)];
axes('Position', cbPos2);
imagesc(clusters_reordered(:));
colormap(gca, clusterColors);
set(gca, 'XTick', [], 'YTick', []);
title('Clusters');

%% Plot indiviudal clusters

nClust = max(cluster_idx);
h = histcounts(cluster_idx,nClust);
[v, idx] = maxk(h,numClusters);

% plot clusters along with medoid 
ncols = ceil(sqrt(numClusters)); 
nrows = ceil(numClusters / ncols);

% 3. Create the tiled layout
figure;
t = tiledlayout(nrows, ncols);

for i = 1:numel(idx)

    nexttile;
idx2 = find(cluster_idx==idx(i));
sequences = {};
for iseq =1:numel(idx2)
    sequences{iseq} = allCurves_z(idx2(iseq),:);
end

clusIdx{i} = idx2;

plot(allCurves_z(cluster_idx==idx(i),:)', 'Color', [.7 .7 .7])
hold on

average = median(vertcat(sequences{:}));
plot(average, 'k', 'LineWidth', 4) 
ax = gca; ax.XTick = 1:7; xlim([0.5 7.5]);
ylim([-2.2 2.2])
defaultAxesProperties(gca, true)
end




