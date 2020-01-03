Ft.mjdg = categorical(Ft.mjdg, 'Ordinal', false);
Ft.task = categorical(Ft.task, 'Ordinal', false);
Ft.sid = categorical(Ft.sid, 'Ordinal', false);
%Ft.mjdg = double(Ft.mjdg); Ft.task = double(Ft.task); Ft.sid = double(Ft.sid);

% not enough variance difference between subjects for the following RE models 
%  qu1(CR):  2(.0269), 3(.0062); 1(1e-17)
%  qu1(PAS): 3(.13); 1(NA), 2(6e-6)
for i = 1:3
    figure
    [P,ANOVATAB,STATS] = anova1(Ft{Ft.task==i & Ft.mjdg==1, 'qu1'}, ...
        Ft{Ft.task==i & Ft.mjdg==1, 'sid'}, 'off');  % kruskalwallis 
    [COMPARISON, MEANS, H, GNAMES] = multcompare(STATS, 'alpha', .05, ...
        'ctype', 'tukey-kramer')
    i, ANOVATAB
    COMPARISON % [GNAMES num2cell(MEANS)]
end

% all effects
[P,T,STATS,TERMS] = anovan(Ft{Ft.task==2, 'qu2'}, Ft{Ft.task==2, {'mjdg' 'ftr1' 'ftr2' 'ps_onset_st' 'ps_offset_st'}}, 'continuous', [2 4 5])


% Whole G anova analysis
G = G7(G7.mjdg==1 & G7.ftr1~=1, :);
fc = {'sid', 'ftr1', 'noi'};
[P,T,STATS,TERMS] = anovan(G.da, G{:, fc}, ...
    'model', 'interaction', 'continuous', 2, 'random', 1, 'varnames', fc)
[rho, pv] = corr(G{:,:}, 'rows', 'pairwise', 'type', 'spearman')
pv < 0.05 / (size(G,2)*(size(G,2)-1)/2)