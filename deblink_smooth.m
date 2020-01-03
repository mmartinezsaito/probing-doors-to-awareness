function [DS, T] = deblink_smooth(P, i, toplt_smoodat)

global fs decf pxPmm

T = double([P{i}{:,1}]); 
D0 = double([P{i}{:,2}]); 
X = T; Y = D0;

%% interpolate blinks
DB = stublinks(D0/100/pxPmm, 0, 0, [], 0.1);
D1 = DB.NoBlinksUnsmoothed;

%{
% reject artifacts by Hampel filtering
D2 = hampel(D1, 3, 3); % only from R2017a?  %  D2 = medfilt1(D1, 3);
[D2, I, Y0, LB, UB, ADX, NO] = hampel_nielsen(T, D1);
figure
plot(T, D1, 'b.'); hold on;      % Original Data
plot(T, D2, 'r');               % Hampel Filtered Data
plot(T, Y0, 'b--');             % Nominal Data
plot(T, LB, 'r--');             % Lower Bounds on Hampel Filter
plot(T, UB, 'r--');             % Upper Bounds on Hampel Filter
plot(T(I), D1(I), 'ks');         % Identified Outliers
%}
D2 = D1;

%% smooth the data with a 25ms Hanning window
wl = round(1000 / fs / decf * 6); % for 240Hz, wl=25ms
DS = conv(D2, hann(wl), 'same')/sum(hann(wl));
%pwelch(hann(wl), [], [], [], fs, 'centered');
%fvtool(hann(25));
if toplt_smoodat  %function pltS(X, Y, M, i_events)
  fh = figure; hold on
  md = median(DS);
  plot(T, DS, 'b-')
  ah = gca;
  axis([get(ah, 'XLim') md-2 md+2])
  line(ones(2,1)*double([M{i}{find(i_events{i} == 1), 1}]), get(ah, 'YLim'), 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
  line(ones(2,1)*double([M{i}{find(i_events{i} == 2), 1}]), get(ah, 'YLim'), 'Color', 'r', 'LineWidth', 1, 'LineStyle', ':');
  plot([M{find(i_events{i}==0), 1}], md, 'ko');
end

