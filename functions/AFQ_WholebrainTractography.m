function fg = AFQ_WholebrainTractography(dt, run_mode, params)
% Perform whole brain deterministic tractography within a white matter mask
%
%      fg = AFQ_WholebrainTractography(dt, [run_mode], params)
%
% Inputs:
% dt       = a dt6 structure (generated by dtiInit).
% run_mode = If 'test' then fewer (10K) fibers are tracked for efficient
%            debugging. This is achieved by setting a higher fa mask
%            threshold and fewer seeds per voxel.  Otherwise don't define.
% params   = A structure of parameters defining the fiber parameters for
%            fiber tracking. Can also be an AFQ structure which contains
%            these parameters
% Outputs:
% fg       = A wholebrain fiber group.
%
% Example:
%
% dt = dtiLoadDt6('subj1/dt6.mat')
% fg = AFQ_WholebrainTractography(dt);
%
% (c) Jason Yeatman, VISTA Lab, Stanford University, 2011

%% Initialize parameters

% We set some reasonable defaults if a parameters structure was not passed
% into the fuction
if ~exist('params','var') || isempty(params)
    
    % Maybe it is better that when there is no defaults, stop it. 
    % Maybe I am not aware that my parameters are not arriving properly,
    % I think it is better to receive an error, and to pass all the parameters
    % always in the config.json files, then you can check what the params
    % are and to try to understand them, otherwise it just works
    %{
    % This defines which voxels will be seeded for tractography.
    if ~exist('run_mode','var') || isempty(run_mode)
        opts.faMaskThresh = 0.30;
    elseif strcmp(run_mode,'test')
        % Set a higher fa threshold, to reduce the number of fibers.
        opts.faMaskThresh = 0.35;
    end
    % Distance between steps in the tractography algoithm
    opts.stepSizeMm = 1;
    % Stopping criteria FA<0.2
    opts.faThresh = 0.2;
    % Discard Fibers shorter than 50mm or longer than 250mm
    opts.lengthThreshMm = [50 250];
    % Stopping criteria angle between steps >30 degrees
    opts.angleThresh = 30;
    % Unknown.....
    opts.wPuncture = 0.2;
    % There are multip.e algorithms that can be used.  We prefer STT. See:
    % Basser PJ, Pajevic S, Pierpaoli C, Duda J, Aldroubi A. 2000.
    % In vivo fiber tractography using DT-MRI data.
    % Magnetic Resonance in Medicine 44(4):625-32.
    opts.whichAlgorithm = 1;
    % Interpolation method. After each step we interpolate the tensor at that
    % point. Trilinear interpolation works well.
    opts.whichInterp = 1;
    % This adds some randomness to each seed point. Each seed point is move
    % randomly by randn*.1mm
    opts.offsetJitter = 0;
    % We seed in voxel in multiple locations. [0.25 and 0.75] Seeds each voxel
    % at 8 equidistant locations.  For one seed in the middle of the voxel use
    % opts.seedVoxelOffsets = 0.5;
    if ~exist('run_mode','var') || isempty(run_mode)
        opts.seedVoxelOffsets = [0.25 0.75];
    elseif strcmp(run_mode,'test')
        opts.seedVoxelOffsets = [0.5];
    end
    % for mrTrix
    % Tracking nfibers and algo to use with mrTrix (the default is for mrTrix3)
    opts.nfibers = 500000;
    opts.mrTrixAlgo = 'iFOD2';
    %}
    error('Whole Brain TRact. is not receiving params, review it please')
elseif isafq(params)
    % If an afq structure is passed in get the tracking parameters and also
    % check if tracking should be done based on CSD with mrtrix
    opts = AFQ_get(params,'tracking parameters');
    mrtrix = AFQ_get(params,'use mrtrix');
    mrtrixVersion = AFQ_get(params,'mrtrixVersion');
else
    % Check to make sure the user defined all the propper parameters
    check = isfield(params,'seedVoxelOffsets')&&isfield(params,'whichInterp')&&...
        isfield(params,'faThresh')&&isfield(params,'whichAlgorithm');
    if check == 0
        error('You did not define all the parameters! Try running with dafaults')
    end
end

%% Track with mrtrix if the right files are there
if exist('mrtrix','var') && mrtrix == 1
    mrtrixpaths = AFQ_get(params,'mrtrix paths',params.currentsub);
    % TODO: make using LiFE an option. Right now, the received fg will
    % have LiFE run on it and the fibers with weight = 0 will be removed.
    [status, results, fg, pathstr] = AFQ_mrtrix_track(mrtrixpaths, ... 
                                                      mrtrixpaths.wm, ... % roi, (wm = wmMask)
                                                      mrtrixpaths.wm_dilated,...  % mask, (wm = wmMask)
                                                      opts.mrTrixAlgo, ...
                                                      opts.nfibers,...
                                                      [],...
                                                      [],...
                                                       1, ...
                                                       mrtrixVersion, ...,
                                                       opts.multishell, ...
                                                       opts.mrtrix_useACT, ...
                                                       opts.life_runLife, ...
                                                       opts.life_discretization, ...
                                                       opts.life_num_iterations, ...
                                                       opts.life_test, ...
                                                       opts.life_saveOutput, ...
                                                       opts.life_writePDB, ...
                                                       opts.ET_runET, ...
                                                       opts.ET_numberFibers, ...
                                                       opts.ET_angleValues, ...
                                                       opts.ET_minlength, ...
                                                       opts.ET_maxlength);
else
    %% Otherwise track with mrdiffusion
    % Compute FA at every voxel
    fa = dtiComputeFA(dt.dt6);
    % Sometimes noise can cause impossible FA values so we clip them
    fa(fa>1) = 1; fa(fa<0) = 0;
    
    %% Create an ROI for wholebrain tractography
    roiAll = dtiNewRoi('all');
    % Make a mask of voxels with FA>faMaskThresh
    mask = fa >= opts.faMaskThresh;
    % Convert mask image to a list of coordinates
    [x,y,z] = ind2sub(size(mask), find(mask));
    clear mask fa;
    % Transofrm mask coordinates to subjects ACPC space
    roiAll.coords = mrAnatXformCoords(dt.xformToAcpc, [x,y,z]);
    clear x y z;
    % Smooth the ROI and fill holes
    roiAll = dtiRoiClean(roiAll,3,{'fillHoles'});
    
    %% Perform wholebrain tractography
    fg = dtiFiberTrack(dt.dt6, roiAll.coords, dt.mmPerVoxel, dt.xformToAcpc, 'wholeBrain', opts);
    
end

return
