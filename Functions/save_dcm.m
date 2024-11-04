function [] = save_dcm(output_folder, d_name, DCM)
    full_path = fullfile(output_folder, [d_name '.mat']);  % Ensure .mat extension
    
    save(full_path, 'DCM');
end
