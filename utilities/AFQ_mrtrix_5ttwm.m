function [status, results] = AFQ_mrtrix_5ttwm(tt5File, ...
                                              faFile, ...
                                              faMaskThresh, ...
                                              brainmask_eroded, ...
                                              wmMask, ...
                                              wmMask_dilated, ...
                                              verbose, ...
                                              mrtrixVersion)

% This is prepared to create the wmmask based on the anatomical 5tt file. 
% The problem is that if we do not have the reversed phase encoding file
% the alligment is not very good so the wmMask we are going to create it
% wont be perfect either...
% 
% Strategy;
% - Use FA > faMaskThresh to create a wmMask for seeding. 
% - Use eroded twice brainmask to clean the borders, it is usually messy
% - Use the wmMask obtained from 5tt file (anatomical) to create the mask
%   where the fibers cannot cross. 
% 
% 
% This is hte code to extract the wm only from the 5tt file
% 
% For now use the FA file to create the mask
%
%  GLU 01.2019


if notDefined('verbose')
    verbose = true;
end
if notDefined('bkgrnd')
    bkgrnd = false;
end


if mrtrixVersion == 2
    error('mrTrix version 2 is deprecated')
end
if mrtrixVersion == 3
    % Create the fa based wmMask first, this will be for seeding
    % Using the eroded brain mask, clean the borders and save with the same name
    cmd_str1 = ['mrthreshold -force -abs ' num2str(faMaskThresh) ' ' ...
                        faFile ' - | mrcalc -force - ' ...
                        brainmask_eroded ' -mult ' wmMask];
        
    % Now create the anatomically based wmMask and dilate it, it will be used to 
    % contain the fibers, but as being it dilated it will encounter the GM
    cmd_str2 = ['mrconvert -coord 3 2 -axes 0,1,2 ' tt5File ' ' ...
                 '- | maskfilter -force -npass 2 - dilate ' ...
                  wmMask_dilated];
    
end

% Send it to mrtrix:
[status,results] = AFQ_mrtrix_cmd(cmd_str1, bkgrnd, verbose,mrtrixVersion);
[status,results] = AFQ_mrtrix_cmd(cmd_str2, bkgrnd, verbose,mrtrixVersion);


