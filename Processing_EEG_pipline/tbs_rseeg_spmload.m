function [spm] = tbs_rseeg_spmload(subject, session, reref)

for i = 1:length(reref)
    n = length(reref(i).time); 
    %m = size(reref(1).time(1));
    
    % Initialize the output matrix
    new_times = {};
    
    % Normalize each epoch
    for j = 1:n
        new_times{j} = linspace(0, 2, 512);
    end
    
    spm(i) = reref(i);
    spm(i).time = new_times;
    spm(i).hdr.chantype(63) = {'EEG'};
    spm(i).hdr.chantype(64) = {'EEG'};
    spm(i).hdr.chantype(65) = {'ECG'};
    spm(i).hdr.chantype(66) = {'EOG'};
    spm(i).hdr.chantype(67) = {'EOG'};
    spm(i).hdr.chantype(68) = {'Marker'};
    spm(i).hdr.chanunit(67) = {'uV'};
    spm(i).hdr.chanunit(68) = {'V'};

    % Generate a unique filename for the SPM file (no extension)
    spm_filename = sprintf('sub-%02d_ses-%02d_epoch-%02d_converted_data', subject, session, i);
    
    % Convert to SPM object and save
    spm_object = spm_eeg_ft2spm(spm(i), spm_filename);
    fprintf('SPM object saved for subject %02d, session %02d, epoch %02d as %s\n', subject, session, i, spm_filename);
end