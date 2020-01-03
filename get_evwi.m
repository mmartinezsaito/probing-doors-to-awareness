function [timelog, rti_qu1, rti_qu2] = get_evwi(hiss, i, s, varargin)

%% extract event window lengths
lkev = {'onset_fp', 'onset_st', 'offset_st', 'onset_qu1', 'rt1', 'onset_qu2', 'rt2'};
if nargin == 4, nz = varargin{1}; end
if s > 6 && nz == 2
  timelog = hiss(i+3).his.timelog; 
else
  timelog = hiss(i).his.timelog;
end
rti_qu1 = timelog(7,:) - timelog(6,:);
rti_qu2 = timelog(9,:) - timelog(8,:); 
