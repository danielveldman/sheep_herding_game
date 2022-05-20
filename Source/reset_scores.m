clear all
close all
clc

% create/clear files for optimization scores
levels = 4;
for ii = levels
    filename = ['optimization_results_level', num2str(ii)];
    Jopt = Inf;
    % TODO check for optimized solutions
    save(filename, 'Jopt');
end

% create/clear files for player scores
levels = 4;
for ii = levels
    filename = ['scores_level', num2str(ii)];
    scores = [];
    optimizer_ind = [];
    % TODO check for optimized solutions
    save(filename, 'scores', 'optimizer_ind');
end
