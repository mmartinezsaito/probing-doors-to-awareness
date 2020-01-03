function [M1, M2, W, ni] = evlk_apd(lkev_str, T, D, M, bls, rws, i_events, ...
                                timelog, facvec, i, fi)
% stimulus onset-locked treatment-averaged pupil diameter data

print_se = 1;
global fs decf

is37 = iscell(T);

if ~is37
  t_onsetst = double([M{i}{find(i_events{i}==0), 1}]);
  t_0 = double([M{i}{find(i_events{i}==1), 1}]);
  t_end = double([M{i}{find(i_events{i}==2), 1}]);
  switch lkev_str
  case 'onset_fp'
    i_onset_fp = timelog(3, :);	
    t_ev = i_onset_fp * 1000 + t_0(2);
  case 'onset_st'
    i_onset_st = timelog(4, :);	
    t_ev = i_onset_st * 1000 + t_0(2);
    %t_ev1 = t_onsetst(find(t_onsetst > t_0(2) & t_onsetst < t_end(2)));
    %norm(t_ev - t_ev1, 2);
  case 'offset_st'
    i_offset_st = timelog(5, :);	
    t_ev = i_offset_st * 1000 + t_0(2);
  case 'onset_qu1'
    i_onset_qu1 = timelog(6, :);	
    t_ev = i_onset_qu1 * 1000 + t_0(2);
  case 'rt1'
    i_rt1 = timelog(7, :);	
    t_ev = i_rt1 * 1000 + t_0(2);
  case 'onset_qu2'
    i_onset_qu2 = timelog(8, :);	
    t_ev = i_onset_qu2 * 1000 + t_0(2);
  case 'rt2'
    i_rt2 = timelog(9, :);	
    t_ev = i_rt2 * 1000 + t_0(2);
  end
  levs = unique(facvec);
  nl = length(levs);
  for li = 1:nl
    ni(li) = sum(facvec == levs(li));
  end
  ne = length(t_ev); median(diff(t_ev));
  assert(ne == length(facvec))
else
  ne = {}; levs = {}; nl = {};
  for j = 1:2
    i2 = i + 3 * (j - 1);
    t_onsetst{j} = double([M{i2}{find(i_events{i2}==0), 1}]);
    t_0{j} = double([M{i2}{find(i_events{i2}==1), 1}]);
    try 
      t_end{j} = double([M{i2}{find(i_events{i2}==2), 1}]);
    catch ME
      warning(ME.identifier, ME.message)
      warning(['The above exception was handled by taking the last '...
               'timepoint of messages log M'])
      tmp = find(i_events{i2}==2); tmp(end) = tmp(end) - 1;
      t_end{j} = double([M{i2}{tmp, 1}]);
    end
    switch lkev_str
    case 'onset_fp'
      i_onset_fp{j} = timelog{j}(3, :);	
      t_ev{j} = i_onset_fp{j} * 1000 + t_0{j}(2);
    case 'onset_st'
      i_onset_st{j} = timelog{j}(4, :);	
      t_ev{j} = i_onset_st{j} * 1000 + t_0{j}(2);
    case 'offset_st'
      i_offset_st{j} = timelog{j}(5, :);	
      t_ev{j} = i_offset_st{j} * 1000 + t_0{j}(2);
    case 'onset_qu1'
      i_onset_qu1{j} = timelog{j}(6, :);	
      t_ev{j} = i_onset_qu1{j} * 1000 + t_0{j}(2);
    case 'rt1'
      i_rt1{j} = timelog{j}(7, :);	
      t_ev{j} = i_rt1{j} * 1000 + t_0{j}(2);
    case 'onset_qu2'
      i_onset_qu2{j} = timelog{j}(8, :);	
      t_ev{j} = i_onset_qu2{j} * 1000 + t_0{j}(2);
    case 'rt2'
      i_rt2{j} = timelog{j}(9, :);	
      t_ev{j} = i_rt2{j} * 1000 + t_0{j}(2);
    end
    ne{j} = length(t_ev{j}); median(diff(t_ev{j}));
    levs{j} = unique(facvec{j});
    nl{j} = length(levs{j});
    for li = 1:nl{j}
      ni(li,j) = sum(facvec{j} == levs{j}(li));
    end
    assert(ne{j} == length(facvec{j})); 
  end
  switch fi
  case 5,    levs = [0 1];  
  case 4,    levs = [-4:-1 1:4];  
  otherwise, levs = union(levs{1}, levs{2}); 
  end 
  nl = length(levs);
end

W = {}; M1 = cell(1,nl); M2 = cell(1,nl);
M1(:) = {0}; M2(:) = {0};
if ~is37
  for k = 1:ne
    o = find(T >= t_ev(k), 1);  % should i have used as baseline t_onset always? it doesn't matter if i only compute differences
    if o - bls < 1 
      blm = mean(D(1:o));
    else
      blm = mean(D(o-bls:o));
    end
    try
      W{k} = (D(o:o+rws-1) - blm) / blm;  % normally DS
    catch ME
      warning(ME.identifier, ME.message)
      warning('The above exception was handled by discarding the window')
      fprintf(1, 'of trial number %d of %d\n', k, ne)
      continue
    end
    f = find(levs == facvec(k));
    M1{f} = M1{f} + W{k};
    M2{f} = M2{f} + W{k}.^2;
  end
else
  for j = 1:2
    for k = 1:ne{j}
      o = find(T{j} >= t_ev{j}(k), 1); 
      blm = mean(D{j}(o-bls:o));
      try           
        W{j,k} = (D{j}(o:o+rws-1) - blm) / blm;  
        if isempty(W{j,k}), error('KuvikException:MissingValue', ...
                'Event %d data is faulty',k), end
      catch ME
        warning(ME.identifier, ME.message)
        warning(['The above exception was handled by zero-padding the ' ...
                 'window of trial number %d of %d\n'], k, ne{j}')
        W{j,k} = zeros(1, rws);
        continue
      end
      f = find(levs == facvec{j}(k));
      M1{f} = M1{f} + W{j,k};
      M2{f} = M2{f} + W{j,k}.^2;
    end
  end 
end

