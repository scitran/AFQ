function [status, results, fg, pathstr] = AFQ_mrtrix_track(files, ...
    roi, ...
    mask, ...
    algo, ...
    nSeeds, ...
    bkgrnd, ...
    verbose, ...
    clobber, ...
    mrtrixVersion, ...
    multishell, ...
    life_runLife, ...    
    life_discretization, ...
    life_num_iterations, ...
    life_test, ...
    life_saveOutput, ...
    life_writePDB, ...
    ET_runET, ...
    ET_numberFibers, ...
    ET_angleValues, ...
    ET_minlength, ...
    ET_maxlength)
%
% function [status, results, fg, pathstr] = mrtrix_track(files, roi, mask, algo, nSeeds, bkgrnd, verbose)
%
% Provided a csd estimate, generate estimates of the fibers starting in roi
% and terminating when they reach the boundary of mask
%
% Parameters
% ----------
% files: structure, containing the filenames generated using mrtrix_init.m
% roi: string, filename for a .mif format file containing the ROI in which
%      to place the seeds. Use the *_wm.mif file for Whole-Brain
%      tractography.
% mask: string, filename for a .mif format of a mask. Use the *_wm.mif file for Whole-Brain
%      tractography.
% algo: Tracking algorithm: it was 'mode' before. Specify it directly in afq.param.track.mrTrixAlgo
% nSeeds: The number of fibers to generate.
% bkgrnd: on unix, whether to perform the operation in another process
% verbose: whether to display standard output to the command window.
% clobber: Whether or not to overwrite the fiber group if it was already
%          computed
%
% Franco, Bob & Ariel (c) Vistalab Stanford University 2013
% Edit GLU 06.2016 added mrTrix versioning
% Mode now has been modified with algo, and it is set in the afq structure
% directly in AFQ_Create. Be careful, the options for mrTrix2 and mrTrix3 are
% very different. See AFQ_Create for all the available options.
% Edit GLU 06.2018 Added Life: after we get the tck we run it through the LiFE pipeline, so that 
% the fiber cleaning is done in origin. Then it is the cleaned WholeBrain
% converted to pdb and continues the normal afq pipeline.
% Edit GLU 08.2018: added Ensemble Tractography

status = 0; results = [];
if notDefined('verbose'),  verbose = false;end
if notDefined('bkgrnd'),    bkgrnd = false;end
if notDefined('clobber'),  clobber = false;end
if notDefined('mrtrixVersion'),    mrtrixVersion = 3;end
if notDefined('multishell'),  multishell = false;end
% LiFE
if notDefined('life_runLife'), life_runLife = true; end
if notDefined('life_discretization'), life_discretization = 360; end
if notDefined('life_num_iterations'), life_num_iterations = 100; end
if notDefined('life_test'), life_test = false; end
if notDefined('life_saveOutput'), life_saveOutput = false; end
if notDefined('life_writePDB'), life_writePDB = false; end
% Ensemble Tractography
if notDefined('ET_runET'), ET_runET = true; end
if notDefined('ET_numberFibers'), ET_numberFibers = 200000; end
if notDefined('ET_angleValues'), ET_angleValues = [47.2, 23.1, 11.5, 5.7, 2.9]; end
if notDefined('ET_minlength'), ET_minlength = 10; end
if notDefined('ET_maxlength'), ET_maxlength = 250; end


% Choose the tracking mode (probabilistic or stream)
% -algorithm name
%      specify the tractography algorithm to use. Valid choices are: FACT, iFOD1,
%      iFOD2, Nulldist1, Nulldist2, SD_Stream, Seedtest, Tensor_Det, Tensor_Prob
%      (default: iFOD2).
% switch mode
%   case {'prob','probabilistic tractography'}
%     mode_str3 = 'Tensor_Prob';
%   case {'stream','deterministic tractogrpahy based on spherical deconvolution'}
%     mode_str3 = 'SD_Stream';
%   case {'tensor','deterministic tractogrpahy based on a tensor model'}
%     mode_str3 = 'Tensor_Det';
%   otherwise
%     error('Input "%s" is not a valid tracking mode', mode);
% end



% This version goes on a gear with mrtrix3 installed, commenting the
% mrtrix2 code and delete it next time
mrtrixVersion = 3;
%{
% In this case there are several changes between mrtrix2 and mrtrix3
% Maintaining the whole thing for version 2 and create a new one for 3 so that
% it will be easier in the future to delete the whole mrTrix2 thing


if mrtrixVersion == 2
    funcName = 'streamtrack';
    % Generate a UNIX command string.
    switch algo
        case {'SD_PROB', 'SD_STREAM'}
            % Build a file name for the tracks that will be generated.
            % THe file name will contain information regarding the files being used to
            % track, mask, csd file etc.
            [~, pathstr] = strip_ext(files.csd);
            tck_file = fullfile(pathstr,strcat(strip_ext(files.csd), '_' , strip_ext(roi), '_',...
                strip_ext(mask) , '_', algo, '-',num2str(nSeeds),'.tck'));
            
            % Generate the mrtrix-unix command
            cmd_str = sprintf('%s %s %s -seed %s -mask %s %s -num %d', ...
                funcName,algo, files.csd, roi, mask, tck_file, nSeeds);
            
        case {'DT_STREAM'}
            % Build a file name for the tracks that will be generated.
            % THe file name will contain information regarding the files being used to
            % track, mask, csd file etc.
            [~, pathstr] = strip_ext(files.dwi);
            tck_file = fullfile(pathstr,strcat(strip_ext(files.dwi), '_' , strip_ext(roi), '_',...
                strip_ext(mask) , '_', algo, '-',num2str(nSeeds),'.tck'));
            
            % Generate the mrtrix-unix command.
            cmd_str = sprintf('%s %s %s -seed %s -grad %s -mask %s %s -num %d', ...
                funcName,algo, files.dwi, roi, files.b, mask, tck_file, nSeeds);
        otherwise
            error('Input "%s" is not a valid tracking algorithm in mrTrix2', algo);
    end
end
%}
% We added ET and LiFE here. 
% The combinations that make sense are:
%    - No ET, No LiFE
%    - No ET, Yes LiFE
%    - Yes ET, Yes LiFE

% Generate a UNIX command string.
switch algo
    case {'iFOD2'}
        % Build a file name for the tracks that will be generated.
        % THe file name will contain information regarding the files being used to
        % track, mask, csd file etc.
        [~, pathstr] = strip_ext(files.csd);
        tck_file = fullfile(pathstr,strcat(strip_ext(files.csd), '_' , strip_ext(roi), '_',...
            strip_ext(mask) , '_', algo, '-',num2str(nSeeds),'.tck'));

        % Generate the mrtrix-unix command
        % See examples at the end of this file
        if  ~multishell
            cmd_str = ['tckgen ' files.csd ' ' ...
                '-algo ' algo ' ' ...
                '-seed_image ' roi ' ' ...
                '-mask ' mask ' ' ...
                '-select ' num2str(nSeeds) ' ' ... % Change on mrtrix interface
                tck_file ' ' ...
                '-force'];
        else
            cmd_str = ['tckgen ' files.csd ' ' ...
                '-algo ' algo ' ' ...
                '-seed_image ' roi ' ' ...
                '-act ' files.tt5 ' ' ...
                '-select ' num2str(nSeeds) ' ' ... % Change on mrtrix interface
                tck_file ' ' ...
                '-force']; % -cutoff 0.05 changed amp thresh cutoff from .1 to 0.05 for testing
        end
    case {'iFOD1'}
        % Build a file name for the tracks that will be generated.
        % THe file name will contain information regarding the files being used to
        % track, mask, csd file etc.
        [~, pathstr] = strip_ext(files.csd);
        tck_file = fullfile(pathstr,strcat(strip_ext(files.csd), '_' , strip_ext(roi), '_',...
            strip_ext(mask) , '_', algo, '-',num2str(nSeeds),'.tck'));

        % Generate the mrtrix-unix command
        % See examples at the end of this file
        if  ~multishell
            cmd_str = ['tckgen ' files.csd ' ' ...
                '-algo ' algo ' ' ...
                '-seed_image ' roi ' ' ...
                '-mask ' mask ' ' ...
                '-num ' num2str(nSeeds) ' ' ...
                tck_file ' ' ...
                '-force'];
        else
            cmd_str = ['tckgen ' files.csd ' ' ...
                '-algo ' algo ' ' ...
                '-seed_image ' roi ' ' ...
                '-act ' files.tt5 ' ' ...
                '-num ' num2str(nSeeds) ' ' ...
                tck_file ' ' ...
                '-force'];
        end
    case {'SD_Stream'}
        % Build a file name for the tracks that will be generated.
        % THe file name will contain information regarding the files being used to
        % track, mask, csd file etc.
        [~, pathstr] = strip_ext(files.csd);
        tck_file = fullfile(pathstr,strcat(strip_ext(files.csd), '_' , strip_ext(roi), '_',...
            strip_ext(mask) , '_', algo, '-',num2str(nSeeds),'.tck'));

        % Generate the mrtrix-unix command
        % See examples at the end of this file
        if  ~multishell
            cmd_str = ['tckgen ' files.csd ' ' ...
                '-algo ' algo ' ' ...
                '-seed_image ' roi ' ' ...
                '-mask ' mask ' ' ...
                '-num ' num2str(nSeeds) ' ' ...
                tck_file ' ' ...
                '-force'];
        else
            cmd_str = ['tckgen ' files.csd ' ' ...
                '-algo ' algo ' ' ...
                '-seed_image ' roi ' ' ...
                '-act ' files.tt5 ' ' ...
                '-num ' num2str(nSeeds) ' ' ...
                tck_file ' ' ...
                '-force'];
        end
    case {'FACT', 'Nulldist1', 'Nulldist2', 'Seedtest', 'Tensor_Det', 'Tensor_Prob'}
        error('Wrapper for algo "%s" not implemented yet. Add the case here and remove it from this list. Use the examples below to build your algorithm.', algo);
    otherwise
        error('Input "%s" is not a valid tracking algo in mrTrix3', algo);
end




if ET_runET && (mrtrixVersion == 3) && strcmp(algo, 'iFOD2') && ~multishell
    disp('Running Ensemble Tractography with mrTrix3, iFOD2, not multishell data.');
    numconcatenate = [];
    for na=1:length(ET_angleValues)
        fgFileName{na}=['fibs' num2str(ET_numberFibers) '_angle' strrep(num2str(ET_angleValues(na)),'.','p') '.tck'];
        fgFileNameWithDir{na}=fullfile(fileparts(tck_file), fgFileName{na});
        cmd_str = ['tckgen ' files.csd ' ' ...
                    '-algo ' algo ' ' ...
                    '-seed_image ' roi ' ' ...
                    '-mask ' mask ' ' ...
                    '-minlength ' num2str(ET_minlength) ' ' ...
                    '-maxlength ' num2str(ET_maxlength) ' ' ...
                    '-angle ' num2str(ET_angleValues(na)) ' ' ...
                    '-select ' num2str(ET_numberFibers) ' ' ...
                    fgFileNameWithDir{na} ' ' ...
                    '-force'];
        % Run it 
        [status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion);
        numconcatenate = [numconcatenate, ET_numberFibers];
    end
    
    % Read, merge, write .mat and .tck file
    % [PATHSTR,NAME,EXT] = fileparts(tck_file);
    % fname = 'NHP_pUM_ETCall_8million_cand.mat';
    et_concatenateconnectomes(fgFileNameWithDir, tck_file, numconcatenate, 'tck'); 
    
else
    fprintf('Running default tracking without Ensemble tractography and mrTrix%d\n', mrtrixVersion);
    
    % Track using the command in the UNIX terminal
    if ~(exist(tck_file,'file') ==2)  || clobber == 1
        [status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose, mrtrixVersion);
    else
        fprintf('\nFound fiber tract file: %s.\n Loading it rather than retracking', tck_file)
    end 
end


% Edit GLU: This is the new part, here we will "clean" the tractogram
% before using it further
% First, obtain the dirNames the mainLife that we already had was using, so
% that we do not change much in the first iteration

mrtrixFolderParts  = split(files.csd, filesep);
% Obtain the session name. This is usually the zip name if it has not
% been edited. 
mrtrixDir  = strjoin(mrtrixFolderParts(1:(length(mrtrixFolderParts)-1)), filesep);
dtiDir     = strjoin(mrtrixFolderParts(1:(length(mrtrixFolderParts)-2)), filesep);
sessionDir = strjoin(mrtrixFolderParts(1:(length(mrtrixFolderParts)-3)), filesep);
lifedir    = fullfile(dtiDir, 'LiFE');

config.dtiinit             = dtiDir;
config.track               = tck_file;
config.life_discretization = life_discretization;
config.num_iterations      = life_num_iterations;
config.test                = life_test;

if life_runLife
    % Change dir to LIFEDIR so that it writes everything there
    if ~exist(lifedir); mkdir(lifedir); end;
    cd(lifedir)

    disp('loading dt6.mat')
    disp(['Looking for file: ' fullfile(config.dtiinit, 'dt6.mat')])
    dt6 = load(fullfile(config.dtiinit, 'dt6.mat'))
    [~,NAME,EXT] = fileparts(dt6.files.alignedDwRaw);
    aligned_dwi = fullfile(sessionDir, [NAME,EXT])

    [ fe, out ] = life(config, aligned_dwi);

    out.stats.input_tracks = length(fe.fg.fibers);
    out.stats.non0_tracks = length(find(fe.life.fit.weights > 0));
    fprintf('number of original tracks	: %d\n', out.stats.input_tracks);
    fprintf('number of non-0 weight tracks	: %d (%f)\n', out.stats.non0_tracks, out.stats.non0_tracks / out.stats.input_tracks*100);

    if life_saveOutput
        disp('writing outputs')
        save('LiFE_fe.mat' ,'fe' , '-v7.3');
        save('LiFE_out.mat','out', '-v7.3');
    else
        disp('User selected not to write LiFE output')
    end

    % This is what we want to pass around
    fg = out.life.fg;


    % Convert the .tck fibers created by mrtrix to mrDiffusion/Quench format (pdb):
    % We will write both, but we want the cleaned ones to be used by vOF or any
    % other downstream code
    pdb_file_LifeFalse = fullfile(pathstr,strcat(strip_ext(tck_file), '_noLiFE.pdb'));
    pdb_file_LifeTrue = fullfile(pathstr,strcat(strip_ext(tck_file), '.pdb'));
    
    if life_writePDB
        mrtrix_tck2pdb(tck_file, pdb_file_LifeFalse);
        mtrExportFibers(fg, pdb_file_LifeTrue, eye(4)); 
    end
else
    % Convert the .tck fibers created by mrtrix to mrDiffusion/Quench format (pdb):
    pdb_file = fullfile(pathstr,strcat(strip_ext(tck_file), '.pdb'));
    fg = mrtrix_tck2pdb(tck_file, pdb_file);
end



end


% tcken options taken from https://github.com/MRtrix3/mrtrix3/blob/master/testing/tests/tckgen
% tckgen dwi.mif   -algo seedtest    -seed_sphere           0,0,4,4  -number 50000 tmp.tck -force && tckmap tmp.tck -template SIFT_phantom/dwi.mif - | testing_diff_data - tckgen/seed_sphere.mif 1000
% tckgen dwi.mif   -algo seedtest    -seed_image            mask.mif -number 3888 tmp.tck -force && tckmap tmp.tck -template SIFT_phantom/dwi.mif - | testing_diff_data - tckgen/SIFT_phantom_seeds.mif 26
% tckgen dwi.mif   -algo seedtest    -seed_random_per_voxel mask.mif 27 tmp.tck -force && tckmap tmp.tck -template SIFT_phantom/dwi.mif - | testing_diff_data - tckgen/SIFT_phantom_seeds.mif 0.5
% tckgen dwi.mif   -algo seedtest    -seed_grid_per_voxel   mask.mif 3 tmp.tck -force && tckmap tmp.tck -template SIFT_phantom/dwi.mif - | testing_diff_data - tckgen/SIFT_phantom_seeds.mif 0.5
% tckgen peaks.mif -algo fact        -seed_rejection        rejection_seed.mif -number 5000 -minlength 4 -mask SIFT_phantom/mask.mif tmp.tck -force && tckmap tmp.tck -template SIFT_phantom/dwi.mif tmp.mif -force && mrstats tmp.mif -mask SIFT_phantom/upper.mif -output mean > tmp1.txt && mrstats tmp.mif -mask SIFT_phantom/lower.mif -output mean > tmp2.txt && testing_diff_matrix tmp1.txt tmp2.txt 30
% tckgen fods.mif  -algo ifod1       -seed_image            mask.mif -act SIFT_phantom/5tt.mif -number 5000 tmp.tck -force && tckmap tmp.tck -template tckgen/act_terminations.mif -ends_only - | mrthreshold - - -abs 0.5 | testing_diff_data - tckgen/act_terminations.mif 0.5
% tckgen dwi.mif   -algo seedtest    -seed_gmwmi            out.mif -act SIFT_phantom/5tt.mif -number 100000 tmp.tck -force && tckmap tmp.tck -template tckgen/gmwmi_seeds.mif - | mrthreshold - - -abs 0.5 | testing_diff_data - tckgen/gmwmi_seeds.mif 0.5
% tckgen peaks.mif -algo fact        -seed_dynamic          fods.mif -mask SIFT_phantom/mask.mif -number 10000 -minlength 4 -initdir 1,0,0 tmp.tck -nthreads 0 -force && tckmap tmp.tck -template SIFT_phantom/dwi.mif tmp.mif -force && mrstats tmp.mif -mask SIFT_phantom/upper.mif -output mean > tmp1.txt && mrstats tmp.mif -mask SIFT_phantom/lower.mif -output mean > tmp2.txt && testing_diff_matrix tmp1.txt tmp2.txt 50
% tckgen fods.mif  -algo ifod2       -seed_image            mask.mif -mask SIFT_phantom/mask.mif -minlength 4 -number 100 tmp.tck -force
% tckgen fods.mif  -algo sd_stream   -seed_image            mask.mif -mask SIFT_phantom/mask.mif -minlength 4 -number 100 -initdir 1,0,0 tmp.tck -force
% tckgen dwi.mif   -algo tensor_prob -seed_image            mask.mif -mask SIFT_phantom/mask.mif -minlength 4 -number 100 tmp.tck -force
% tckgen fods.mif  -algo ifod1       -seed_image            mask.mif -act SIFT_phantom/5tt.mif -backtrack -number 100 tmp.tck -force
% tckgen dwi.mif   -algo tensor_det  -seed_grid_per_voxel   mask.mif 3 -nthread 0 tmp.tck -force && testing_diff_tck tmp.tck tckgen/tensor_det.tck 1e-2
% tckgen dwi.mif   -algo tensor_det  -seed_grid_per_voxel   mask.mif 3 tmp.tck -force && testing_diff_tck tmp.tck tckgen/tensor_det.tck 1e-2
