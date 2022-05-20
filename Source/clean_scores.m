clear all
close all
clc

levels = 1:4;

for level = levels
    load(['scores_level', num2str(level)])
    
    clear Js
    for jj = 1:length(scores)
        if jj ~= optimizer_ind
            if isequal(scores(jj).name, 'Name')
                Js(jj) = Inf;
            else
                Js(jj) = scores(jj).J;
            end
        end
    end
    
    [Js,ind] = sort(Js);
    scores_old = scores;
    clear scores
    scores(1).J    = scores_old(optimizer_ind).J;
    scores(1).name = scores_old(optimizer_ind).name;
    scores(1).u    = scores_old(optimizer_ind).u;
    optimizer_ind = 1;
    
    for jj = 2:3
        scores(jj).J    = scores_old(ind(jj)).J;
        scores(jj).name = scores_old(ind(jj)).name;
        scores(jj).u    = scores_old(ind(jj)).u;
    end
    
    save(['scores_level', num2str(level)], 'optimizer_ind', 'scores')
end