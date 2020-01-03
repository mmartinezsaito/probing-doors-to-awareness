function [iei, avg_iei, std_iei] = plt_evwist(timelog, toplt)
% event window length averages: fpdi;(2) sdi; onset_qu1 - offset_st;(4) rti1; onset_qu2 - rt1;(6) rti2

rtil = {'fpdi', 'sdi', 'qu1on-stoff', 'rti1', 'qu2on-rt1', 'rti2', 'isi'}; 
iei = diff(timelog(3:9,:) - repmat(timelog(3,:), [7 1])); 
iei(7,:) = [timelog(3,2:end)-timelog(9,1:end-1) timelog(3,2)-timelog(9,1)];
avg_iei = mean(diff(timelog(3:9,:) - repmat(timelog(3,:), [7 1])), 2);  %diff(mean(timelog(3:9,:) - repmat(timelog(3,:), [7 1]), 2));
avg_iei(7) = mean(timelog(3, 2:end) - timelog(9, 1:end-1)); % average: onset_fp(n) - rt2(n-1) 
std_iei = std(diff(timelog(3:9,:) - repmat(timelog(3,:), [7 1]), 1), 0, 2);
std_iei(7) = std(timelog(3, 2:end) - timelog(9, 1:end-1));
if toplt
  figure, boxplot(iei', 'jitter', 0, 'notch', 'on', 'datalim', [0 3], ...
      'extrememode', 'compress', 'labels', rtil)
end
