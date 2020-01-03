


lme1 = 'right1  ~ (ftr1 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2) + (1+ftr1 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2|sid)';
lme1 = 'qu1  ~ (ftr1 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2) + (1+ftr1 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2|sid)';
lme2 = 'qu2  ~ (ftr1 + qu1 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2) + (1+ftr1+qu1+ps_onset_st+ps_onset_qu2+ps_rt2|sid)';
lme2_m3 = 'qu2  ~ (ftr1 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2) + (1+ftr1+ps_onset_st+ps_onset_qu2+ps_rt2|sid)';
lme3 = 'qu12 ~ (ftr1 + rti1 + rti2 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2) + (1|sid)';
lme4 = 'rti1 ~ (ftr1 + qu1 + qu2 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2) + (1+ftr1+qu1+qu2+ps_onset_st+ps_onset_qu1+ps_onset_qu2+ps_rt2|sid)';
lme4a = 'rti1 ~ (ftr1 + right1 + right2 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2) + (1+ftr1+right1+right2+ps_onset_st+ps_onset_qu1+ps_onset_qu2+ps_rt2|sid)';
lme4lr = 'rti1 ~ (ftr1 + tk2qu1lr + qu2 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2) + (1+ftr1+tk2qu1lr+qu2+ps_onset_st+ps_onset_qu1+ps_onset_qu2+ps_rt2|sid)';
lme5 = 'rti2 ~ (ftr1 + rti1 + qu1 + qu2 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2) + (1+ftr1+rti1+qu1+qu2+ps_onset_st+ps_onset_qu1+ps_onset_qu2+ps_rt2|sid)';
lme5a = 'rti2 ~ (ftr1 + rti1 + right1 + right2 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2) + (1+ftr1+rti1+right1+right2+ps_onset_st+ps_onset_qu1+ps_onset_qu2+ps_rt2|sid)';
lme5lr = 'rti2 ~ (ftr1 +rti1 + tk2qu1lr + qu2 + ps_onset_st + ps_onset_qu1 + ps_onset_qu2 + ps_rt2) + (1+ftr1+rti1+tk2qu1lr+qu2+ps_onset_st+ps_onset_qu1+ps_onset_qu2+ps_rt2|sid)';
lme6 = 'rti1 ~ (ftr1        + right1 + right2 + ps_onset_st + ps_onset_qu1 + ps_rt1 + ps_onset_qu2 + ps_rt2) + (1|sid)';
lme7 = 'rti2 ~ (ftr1 +rti1  + right1 + right2 + ps_onset_st + ps_onset_qu1 + ps_rt1 + ps_onset_qu2 + ps_rt2) + (1|sid)';
lme8 = 'right1 ~ (ftr1 + rti1 + rti2 + ps_onset_st + ps_rt1 + ps_onset_qu2 + ps_rt2) + (1 | sid)';
lme9 = 'right2 ~ (ftr1 + rti1 + rti2 + ps_onset_st + ps_rt1 + ps_onset_qu2 + ps_rt2) + (1 | sid)';
lme10 = 'ps_onset_st  ~ (ftr1) + (1|sid)';
lme11 = 'ps_onset_qu1  ~ (ftr1) + (ftr1|sid)';
lme12 = 'ps_rt1  ~ (ftr1 + qu1) + (ftr1+qu1|sid)';
lme13 = 'ps_onset_qu2  ~ (ftr1 + qu1) + (ftr1+qu1|sid)';
lme14 = 'ps_rt2  ~ (ftr1 + qu1 + qu2) + (ftr1+qu1+qu2|sid)';
lme14_m3 = 'ps_rt2  ~ (ftr1 + qu2) + (ftr1+qu2|sid)';

for i = [1 2 3], for m = [1 3],m,i
    if m == 3 && i == 1 &&  1
      lme = lme2_m3; % singular matrix for i=1,m=3 because all(qu1==1)! 
    end
        
    % CR, VR
    %GLME{m, i} = fitglme(Ft(Ft.task==i & Ft.mjdg==m, :), lme2, ...
    %    'Distribution', 'binomial', 'BinomialSize', 4, 'Link', 'logit', ...
    %    'CovariancePattern', 'FullCholesky', 'FitMethod', 'Laplace', ...
    %    'DummyVarCoding', 'reference', 'Optimizer', 'quasinewton', ...
    %    'StartMethod', 'default', 'Verbose', true, 'CheckHessian', true);
    
    % reaction times
    %GLME{m, i} = fitglme(Ft(Ft.task==i & Ft.mjdg==m & ~rt4s_out, :), lme1, ...
    %    'Distribution', 'Gamma', 'Link', 'log', ...
    %    'CovariancePattern', 'FullCholesky', 'FitMethod', 'Laplace', ...
    %    'DummyVarCoding', 'reference', 'Optimizer', 'quasinewton', ...
    %    'StartMethod', 'default', 'Verbose', true, 'CheckHessian', true);
    %LME{m,i} = fitlme(Ft(Ft.task==i & Ft.mjdg==m & ~rt4s_out, :), lme7, ...
    %    'CovariancePattern', 'FullCholesky', 'FitMethod', 'ML', ...
    %    'DummyVarCoding', 'reference', 'Optimizer', 'quasinewton', ...
    %    'StartMethod', 'default', 'Verbose', true, 'CheckHessian', true)
    %LM{m,i} = fitlm(Ft(Ft.task==i & Ft.mjdg==m & ~rt4s_out, :), lme5)                
    %GLM{m,i} = fitglm(Ft(Ft.task==i & Ft.mjdg==m & ~rt4s_out, :), lme4, ...
    %    'Distribution', 'normal', 'Link', 'log')
    
    % pupil amplitudes
    LME{m,i} = fitlme(Ft(Ft.task==i & Ft.mjdg==m, :), lme10, ...
        'CovariancePattern', 'FullCholesky', 'FitMethod', 'ML', ...
        'DummyVarCoding', 'reference', 'Optimizer', 'quasinewton', ...
        'StartMethod', 'default', 'Verbose', true, 'CheckHessian', true)
end, end

fit = GLME{1,1};
plotResiduals(fit)
scatter(fit.fitted, fit.residuals)
coefCI(fit)
randomEffects(fit)
%anova(fit, 'DFMethod', 'satterthwaite')
fit.Rsquared
fit.CoefficientCovariance
%rank(fit.designMatrix)


% meta_da-da dip for detec: a trend p=0.06
G=G4f13(G4f13.task==1&G.mjdg==1,:);
[P,T,STATS,TERMS] = anovan(G.meta_da-G.da, G{:,{'sid' 'ftr1'}},'random',1,'display','off');
[P,ANOVATAB,STATS]=anova1(G.meta_da-G.da, G.ftr1)
[COMPARISON, MEANS, H, GNAMES] = multcompare(STATS, 'alpha', .05, 'ctype', 'tukey-kramer')