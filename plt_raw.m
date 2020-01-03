function plt_raw(P, M, i_events, s, I)

global fs decf 
if nargin < 5
  switch s < 7
  case 1, I = 1:3;
  case 0, I = 1:6;
  end 
end


for i = I 
  fh = figure;
  ah(1) = subplot(2,2,1); hold on
  plot([P{i}{:,1}], [P{i}{:,2}]/100, 'b-')
  ylabel('left pupil diameter (px)')
  ah(2) = subplot(2,2,2); hold on
  plot([P{i}{:,1}], [P{i}{:,3}]/100, 'b-')
  ylabel('left pupil area (px^2)')
  ah(3) = subplot(2,2,3); hold on
  plot([P{i}{:,1}], [P{i}{:,4}]/100, 'm-', [P{i}{:,1}], [P{i}{:,5}]/100, 'g-')
  ylabel('left pupil position (px)'), legend('x', 'y')
  ah(4) = subplot(2,2,4);
  plot([P{i}{:,4}]/100, [P{i}{:,5}]/100, 'c-')
  xlabel('left pupil x position (px)'), ylabel('left pupil y position (px)')
  for l = 1:3
    axes(ah(l))
    title(sprintf('Samplig period: %0.2fms', decf/fs*1000)), xlabel('time (ms)') 
    line(ones(2,1)*double([M{i}{find(i_events{i} == 1), 1}]), get(ah(l), 'YLim'), 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
    line(ones(2,1)*double([M{i}{find(i_events{i} == 2), 1}]), get(ah(l), 'YLim'), 'Color', 'r', 'LineWidth', 1, 'LineStyle', ':');
    plot([M{i}{find(i_events{i}==0), 1}], median([P{i}{:,l+1}])/100, 'ko');
  end
end

