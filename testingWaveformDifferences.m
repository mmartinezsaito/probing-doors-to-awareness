%  printing sample-by-sample massive t-testing
%  fit loads of regression models and anovas 
%  the proper way would be cluster-based permutation test

s = 1, m = 1, i = 1

anatype = 'fitlm'
%regstats
Fm_a = Fm(:, 1:4);
Fm_r = [ones(size(Fm, 1) ,1) Fm_a];

apv = {}; csg = {}; B = {}; Bci = {}; B1 = {}; LM = {};
for ei = 1:length(lkev)
    csg{ei} = 0;
    for t = 1:size(W_smie{s,m,i,ei}, 2)
        switch anatype
        case 'anovan'
            [apv{ei, t}, atab, astats, aterms] = anovan(W_smie{s,m,i,ei}(:,t), Fm_a, ...
                'varnames', dsf(1:4), 'model', 'linear', 'display', 'off');
            disp(apv{ei,t}')
            csg{ei} = csg{ei} + (apv{ei,t} < 0.05 / size(W_smie{s,m,i,ei},2));
        case 'linreg'
            [B{ei,t}, Bci{ei,t}, res, outint, RsqFPv] = regress(W_smie{s,m,i,ei}(:, t), Fm_r);
        case 'ridge'
            B1{ei,t} = ridge(W_smie{s,m,i,ei}(:,t), Fm_a, 0:.1:1);
        case 'fitlm'
            LM{ei,t} = fitlm(Fm_a, W_smie{s,m,i,ei}(:,t), sprintf('PD ~ %s*%s*%s*%s', ...
                dsf{1:4}), 'VarNames', ['PD' dsf(1:4)]);
        end
    end
end

switch anatype
    case 'anovan'
        [csg{:}]'
        % Holm-Bonferroni
        [corrected_p, h] = bonf_holm(cell2mat(apv), 0.05);
        sum(h(:)), length(h(:))
        % FDR
        [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(cell2mat(apv), 0.05, 'pdep', 'yes');
    case 'ridge'
        plot(k, B{1}');
        xlabel('Ridge parameter'); ylabel('Standardized coef.');
        title('Ridge Trace')
        lgdc = arrayfun(@num2str, dsf(1:4), 'UniformOutput', 0);
        legend(lgdc, 'Location', 'Best');
    case 'fitlm'
        LM{1}.plot
end