%────────────────────────────────────────────────
% automatically find all “N_output.dat” files in output/
files = dir('output/*_output.dat');

if isempty(files)
  error('No *_output.dat files found in output/ folder.');
end

% extract the leading number from each filename
nums = arrayfun(@(f) sscanf(f.name, '%d_output.dat'), files);

% guard against any non-matches
nums = nums(~isnan(nums));

% set low and high indexes
low_index  = min(nums);
high_index = max(nums);

% print out all the indices we found
fprintf('Found output files with indices: ');
fprintf('%d ', nums);
fprintf('\n');

fprintf('→ low_index = %d, high_index = %d\n', low_index, high_index);

% now build your range
directory_range = low_index : high_index;
%────────────────────────────────────────────────
