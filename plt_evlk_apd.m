function plt_evlk_apd(AV, SE, rwt, facvec, lkev_str, gbf, print_se)  

levs = unique(facvec);
nl = length(levs);
cm = hsv(nl);
lgdc = arrayfun(@num2str, levs, 'UniformOutput', 0); 
%lgd = cell(1, 3*lf); lgd(:) = {''}; [lgd{1:3:3*nl}] = deal(lgdc{:});
figure, hold on
for li = 1:nl
  ah(li) = plot(rwt, AV{li}, '-');
  set(ah(li), 'Color', cm(li,:)) %for c = ahc, set(c, 'Color', cm(i, :)), end 
  if print_se
    ahc = plot(rwt, AV{li}-SE{li}, '--', rwt, AV{li}+SE{li}, '--');
    for c = ahc, set(c, 'Color', cm(li, :)), end 
  end
end
title(['Grouped by ' gbf ' (' num2str(nl) ' levels: ' sprintf('%d ', levs) ')'])
xlabel(['Time (ms); zero locked to ' strrep(lkev_str, '_', '\_')])
ylabel('Normalized pupil diameter')
legend(ah, lgdc, 'Location', 'Best')
