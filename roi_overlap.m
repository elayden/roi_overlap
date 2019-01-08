function [output_table] = roi_overlap(clusters, atlas, labels)
%AAL_COMPARE Outputs table showing percentage overlap of cluster w/ brain 
%     atlases 
% 
% Note: to work properly, the clusters or cluster file must be in the same
% space as the desired atlas (same dimensions & registered)
% 
% INPUTS:
%   clusters,            1. a cell array with linear clusters of a cluster in each
%                       cell, or x,y,z,subscripts; 2. an image file (.img/.hdr, 
%                       .nii, .nii.gz) wherein clusters are >0 voxel values, 
%                       and all else == 0.
%   atlas,              Integer: 1 = AAL, 2 = Harvard-Oxford, 3 = Brodmann Areas
%                       OR, filepath string specifying custom ROI/atlas file
%                       (.img,.nii,.nii.gz)
%   labels (optional),  use this input only if a custom atlas file is used 
%                       for previous input; should be full filepath string
%                       denoting a .txt file; text file should contain a
%                       list of ROIs for the custom atlas separated by
%                       returns/new lines
% 
% OUTPUTS:
%   output_table,   a cell array table showing labels and percent
%                   overlap
%   indrois,        clusters within atlas
%   perc_reg,       percent of region covered by cluster
%   numvoxels,      number of voxels within region

% Identify Function Path and Add Helper Scripts:
script_fullpath = mfilename('fullpath');
[script_path,~,~] = fileparts(script_fullpath);
addpath(genpath(script_path))

% Load appropriate labels and atlas:
if isnumeric(atlas)
    switch atlas
        case 1,
        label_img = load_nii('aal.nii');
        labels = importdata('aal.txt');
        labels{length(labels)+1} = 'Not Labelled';
        template_name = 'AAL Region';
        case 2,
        label_img = load_nii('HarvardOxford.nii'); 
        labels = importdata('HarvardOxford.txt');
        labels{length(labels)+1} = 'Not Labelled'; 
        template_name = 'Harvard-Oxford Region';
        case 3,
        label_img = load_nii('BA.img'); 
        labels = importdata('BA.txt');
        labels{length(labels)+1} = 'Not Labelled';    
        template_name = 'Brodmann Area';
    end
elseif ischar(atlas)
    try
        label_img = load_nii(atlas); 
    catch
        error('Could not load specified custom atlas file.')
    end
    try
        labels = importdata(labels);
    catch
        error('Could not load specified labels for custom atlas file.')
    end
    labels{length(labels)+1} = 'Not Labelled'; 
    template_name = 'Custom Atlas';
else
    error('Input ''atlas'' must be specified as an integer denoting atlas type or as a string path to a custom atlas image file.')
end
atlas_dim = size(label_img.img);

% Check input type for 'clusters':
if iscell(clusters)
    nClust = length(clusters);
    if size(clusters{1},2)==3 % subscripts
        % Convert to linear indices:
        for i = 1:nClust
            clusters{i} = sub2ind(atlas_dim,clusters{i}(:,1),clusters{i}(:,2),clusters{i}(:,3));
        end
    elseif size(clusters{1},2)~=1
        error('Input ''clusters'', if cell array, should contain either linear indices or 3 dimensional subscripts.')
    end
elseif ischar(clusters)
    try
        clust_img = load_nii(clusters); 
        % Check dimensions:
        if ~all(atlas_dim==size(clust_img.img))
            error('''clusters'' image dimensions must match ''atlas'' image dimensions.')
        end
    catch
        error('Could not load specified ''clusters'' image file.')
    end
    nClust = max(clust_img.img(:));
    % Extract linear indices from 3D clusters image:
    clusters = cell(1,nClust);
    for i = 1:nClust
        clusters{i} = ind2sub(atlas_dim,find(clust_img.img(:)==i));
    end
end

% Create output table:
output_table = cell(1,nClust);
for i = 1:nClust
    rois = label_img.img(clusters{i});
    indrois = unique(rois); % how many atlas ROIs overlap?
    nOverlaps = length(indrois);
    numvoxels = zeros(1,nOverlaps);
    perc_reg = zeros(1,nOverlaps);
    for j = 1:nOverlaps
        numvoxels(j) = sum(rois==indrois(j));
        perc_reg(j) = numvoxels(j)/sum(label_img.img(:)==indrois(j));
    end
    perc_reg = round(perc_reg.*100,2);
    if any(indrois==0); indrois(indrois==0) = length(labels); end
    [numvoxels,ix] = sort(numvoxels,'descend'); 
    indrois = indrois(ix); 
    perc_reg = perc_reg(ix); 
    output_table{i}(1,1:3) = {[template_name,' - ',num2str(i)],'# Voxels','% Region Covered'};
    output_table{i}(2:length(ix)+1,1) = labels(indrois);
    output_table{i}(2:length(ix)+1,2) = num2cell(numvoxels);
    count = 0;
    for j = 2:length(ix)+1
        count = count + 1;
        output_table{i}{j,3} = sprintf('%.3g',perc_reg(count));
    end
end

end % End roi_overlap.m
