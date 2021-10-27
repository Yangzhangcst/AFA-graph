clc;clear ; close all;
%% set path
addpath 'others'
addpath 'evals'
addpath 'OLRSC'
addpath 'SSC'
addpath 'functions'
addpath 'msseg'
addpath 'APC'
addpath 'Graph_based_segment'

%% set parameters for bipartite graph
para.alpha = 0.001; % affinity between pixels and superpixels
para.beta  =  20;   % scale factor in superpixel affinity
para.nb = 1; % number of neighbors for superpixels
para.rho = 1;
para.epochs = 20;
para.d = 50;
para.Nimgs = 591; % number of images in BSDS300/500

%% read numbers of segments used in the paper
bsdsRoot = './database/MSRC/';
saveRoot = 'results_searchN_denoise_n';
fid = fopen(fullfile('Nsegs_MSRC.txt'),'r');
line = 1;
while feof(fid) == 0    
    BSDS_INFO{line,1} = deblank(fgetl(fid)); 
    line = line+1;
end
fclose(fid);

%% PRI,VoI,GCE,BDE.
PRI_all = zeros(para.Nimgs,1);
VoI_all = zeros(para.Nimgs,1);
GCE_all = zeros(para.Nimgs,1);
BDE_all = zeros(para.Nimgs,1);


%% Settings
% setup
nGCluster = 3; % number of subjects.
% dimension reduction
reduceDimension = @(data) dimReduction_PCA(data, 0);
% normalization
normalizeColumn = @(data) cnormalize_inplace(data);
% representation
buildRepresentation = @(data) OMP_mat_func(data, 3, 1e-6); % second parameter is sparsity
% spectral clustering   
genLabel = @(affinity, nCluster) SpectralClustering(affinity, nCluster, 'Eig_Solver', 'eigs');

%%
Nseg_save = [];
for idxI =1:para.Nimgs
    % read number of segments
    tic; img_name = BSDS_INFO{idxI};
    img_loc = fullfile(bsdsRoot,'images',[img_name,'.bmp']);
    present = 'image';
    img = im2double(imread(img_loc)); [X,Y,~] = size(img);
    out_path = fullfile(saveRoot,'MSRC','ikde_1_1_APC_3_5',img_name);
    if ~exist(out_path,'dir'), mkdir(out_path); end
    if 1
    % generate superpixels
    [para_MS, para_FH] = set_parameters_oversegmentation(img_loc);
    [seg,labels_img,seg_vals,seg_lab_vals,seg_edges,seg_img] = make_superpixels(img_loc,para_MS,para_FH);

    %% construct graph
    Np = X*Y;   Nsp = 0;
    for k = 1:length(seg),  Nsp = Nsp + size(seg{k},2); end
    W_Y = sparse(Nsp,Nsp);  edgesXY = [];   j = 1;
    for k = 1:length(seg) 
        % for each over-segmentation
        feature = seg_lab_vals{k};
        feature = ikde(feature,1,1);

        tmp1 = reduceDimension(feature);
        tmp1 = normalizeColumn(tmp1);
        R = buildRepresentation(tmp1');
        R(1:length(feature)+1:end) = 0;
        A = abs(R) + abs(R)';
        nGCluster(idxI,k) = APclustering_w_o_ikde(feature);
        index_tmp = genLabel(A, nGCluster(idxI,k));

        % superpixel division
        local_nodes  = find(index_tmp == mode(index_tmp));
        global_nodes = find(index_tmp ~= mode(index_tmp));

        % first we construct the adjacent graph over all nodes
        w = makeweights(seg_edges{k},feature,para.beta);
        W_local = adjacency(seg_edges{k},w);

        % assign local graph entries to fused new graph W, we will 
        % replace the nodes belongs to globla_nodes with value of
        % global graph value
        %W=zeros(size(feature,1),size(feature,1));
        W=W_local;

        % randomly generate two set of supperpxiels from Medium set
        p = global_nodes;
%         p = randperm(length(global_nodes));
        % please choose different kinds of global graph combination
        W_OLRSC = OLRSCGRAPH_n(feature,para); 
        W = assignGraphValue(W,W_OLRSC,p);
        W = sparse(W);

        Nk = size(seg{k},2); % number of superpixels in over-segmentation k
        W_Y(j:j+Nk-1,j:j+Nk-1) = prune_knn(W,para.nb);

        % affinities between pixels and superpixels
        for i = 1:Nk
            idxp = seg{k}{i}; % pixel indices in superpixel i
            Nki = length(idxp);
            idxsp = j + zeros(Nki,1);
            edgesXY = [edgesXY; [idxp, idxsp]];
            j = j + 1;
        end
    end
    W_XY = sparse(edgesXY(:,1),edgesXY(:,2),para.alpha,Np,Nsp);
    % affinity between a superpixel and itself is set to be the maximum 1.
    W_Y(1:Nsp+1:end) = 1;  B = [W_XY;W_Y];

    %% Transfer cut
    out_path_gt= fullfile(saveRoot,'MSRC','ikde_1_1_APC_3_5',img_name);
    if ~exist(out_path_gt,'dir'), mkdir(out_path_gt); end
    [gt_imgs, gt_cnt] = view_gt_segmentation(bsdsRoot,img,present,out_path_gt,img_name,para,0);
    nclusters=1:40;  E=zeros(length(nclusters),5);segs=cell(1,length(nclusters));
    for ncluster=1:length(nclusters)
        if size(B,2) < nclusters(ncluster), break; end
        label_img = Tcut(B,nclusters(ncluster),[X,Y]);
        % display the result
        view_segmentation(img,label_img(:),out_path,img_name,0);
        % Evaluation and save result
        out_vals = eval_segmentation(label_img,gt_imgs);
        E(ncluster,:)=[nclusters(ncluster),out_vals.PRI,out_vals.VoI,out_vals.GCE,out_vals.BDE];
        segs{ncluster}=uint16(label_img);
    end
    out_seg_path = fullfile(saveRoot,'MSRC','ikde_1_1_APC_3_5','segs');
    if ~exist(out_seg_path,'dir'), mkdir(out_seg_path); end
    out_seg = fullfile(out_seg_path,[img_name, '.mat']);
    save('-v7',out_seg, 'segs');
    outname = fullfile(out_path,[img_name, '.mat']);
    save('-v7',outname, 'E');else
    outname = fullfile(out_path,[img_name, '.mat']);load(outname);end
    % Evaluation and save result
    [maxE,idx] = max(E(:,2));ti = toc;
    fprintf('%6d %6s: %2d %9.6f, %9.6f, %9.6f, %9.6f %.4fs\n', idxI,...
        img_name,E(idx,1), E(idx,2), E(idx,3), E(idx,4), E(idx,5),ti);
    Neg_all(idxI) = E(idx,1);
    PRI_all(idxI) = E(idx,2);  VoI_all(idxI) = E(idx,3);
    GCE_all(idxI) = E(idx,4);  BDE_all(idxI) = E(idx,5);
end

fprintf('Mean: %14.6f, %9.6f, %9.6f, %9.6f \n', mean(PRI_all), mean(VoI_all), mean(GCE_all), mean(BDE_all));
fid_out = fopen(fullfile(saveRoot,'MSRC','ikde_1_1_APC_3_5','evaluation.txt'),'w');
for idxI=1:para.Nimgs
    fprintf(fid_out,'%6s %9.6f, %9.6f, %9.6f, %9.6f \n', BSDS_INFO{idxI},...
        PRI_all(idxI), VoI_all(idxI), GCE_all(idxI), BDE_all(idxI));
end
fprintf(fid_out,'Mean: %10.6f, %9.6f, %9.6f, %9.6f \n', mean(PRI_all), mean(VoI_all), mean(GCE_all), mean(BDE_all));
fclose(fid_out);

fid_out2 = fopen(fullfile(saveRoot,'MSRC','ikde_1_1_APC_3_5','Nsegs.txt'),'w');
for idxI=1:para.Nimgs
    fprintf(fid_out2,'%6s %d \n', BSDS_INFO{idxI},Neg_all(idxI));
end
fclose(fid_out2);
