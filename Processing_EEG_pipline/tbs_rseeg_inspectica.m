function comps = tbs_rseeg_inspectica(comps)
% Identify which components are to be removed
% Done through automated process with follow-up visual confirmation
% Automated process to suggest components for removal
% Plot, inspect, and identify first 25 ICs for removal (Display 5 x 5)
% At the prompt, list the components to reject
% Example: [3,6,7,8,9];

%% Automated ICA inspection
[blinks, eye_movements, muscle, gen_disc, suggested_comps,all_art] = tbs_rseeg_autoica(comps);

%% Visual ICA inspection
tbs_rseeg_plotica(comps,suggested_comps);

%% List components for rejection
% Suggestion of components to remove
disp('Automated ICA reject components are:')
disp(['Blinks: ', num2str(blinks)]);
disp(['Lateral eye movements: ', num2str(eye_movements)]);
disp(['Generic discontinuities: ', num2str(gen_disc)]);
disp(['Muscle: ', num2str(muscle)]);
disp(['All artifacts: ', num2str(all_art)]);
disp('.');
disp(['Suggested: ', num2str(suggested_comps)]);

%% Create prompt for user to include which components need to be rejected
prompt          = 'list components to reject. Approximately 10-15% of total:';
x               = input(prompt);
comps.rejected  = x;
close
end
