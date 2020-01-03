%readPupilData1
%readPupilData2

% load data
pars = {}; hiss = {}; PA = {}; MA = {}; i_eventsA = {}; 
DS = {}; T = {}; DS1 = {}; T1 = {}; 
timelog = {}; timelog1 = {}; rti_qu1 = {}; rti_qu2 = {};
loadDataFilenames

for s = ss 
  for m = ms 
    if     s <=6 
      load([datapath 'E3_4/s' sprintf('%02d', s) 'm' sprintf('%1d', m) 'c']) 
      for i = is
        pars{s,m}(i) = load([datapath 'E3_4/' sfn{s, m + 1}{i}], 'par');
        hiss{s,m}(i) = load([datapath 'E3_4/' sfn{s, m + 1}{i}], 'his');
      end
    elseif s > 6
      load([datapath 'E3_7/s' sprintf('%02d', s) 'm' sprintf('%1d', m) 'c']) 
      for i = is
        for j = [0 3]
          pars{s,m}(i+j) = load([datapath 'E3_7/' sfn{s, m + 1}{i+j}], 'par');
          hiss{s,m}(i+j) = load([datapath 'E3_7/' sfn{s, m + 1}{i+j}], 'his');
        end
      end 
    end % E3_4 or E3_7     
    PA{s,m} = P; % pupil size 
    MA{s,m} = M; i_eventsA{s,m} = i_events;  % time and event logs 
    
    for i = is
      % deblink and smooth
      [DS{s,m,i}, T{s,m,i}] = deblink_smooth(PA{s,m}, i, toplt_smoodat);
      if s > 6
        [DS1{s,m,i}, T1{s,m,i}] = deblink_smooth(PA{s,m}, i+3, toplt_smoodat); 
      end

      % deconvolve the whole DS timeseries? not a good idea
      %DP{s,m,i} = deconvPsts(DS{s,m,i}, rwt * 1000, 0);               
      
      
      % extract event window lengths
      if s <= 6
        [timelog{s,m,i}, rti_qu1{s,m,i}, rti_qu2{s,m,i}] = get_evwi(hiss{s,m},i,s);
      elseif s >= 7
        for nz = 1:2 
          [timelog{s,m,i,nz}, rti_qu1{s,m,i,nz}, rti_qu2{s,m,i,nz}] = get_evwi(hiss{s,m},i,s,nz);
        end  
      end
      
      % divide the signal (0 to 2500~3000ms) by the median pupil size during a baseline period (-200~100 to 0ms)
      if s < 7
        % ftr1: luminance8 or angle4 or luminance3;  ftr2: a_rps3
        qu12 = ((2 * hiss{s,m}(i).his.qu1 - 1) .* hiss{s,m}(i).his.qu2); 
        Fm{s,m,i} = [hiss{s,m}(i).his.ftr(1,:)' hiss{s,m}(i).his.ftr(end,:)' hiss{s,m}(i).his.qu1' hiss{s,m}(i).his.qu2' qu12'];
        Ft_smi{s,m,i} = array2table(Fm{s,m,i}, 'VariableNames', dsf);
      else
        nt0 = pars{s,m}(i).par.nt; nt1 = pars{s,m}(i+3).par.nt;
        qu12 = ((2 * hiss{s,m}(i).his.qu1 - 1) .* hiss{s,m}(i).his.qu2); 
        qu12n1 = ((2 * hiss{s,m}(i+3).his.qu1 - 1) .* hiss{s,m}(i+3).his.qu2);
        
        % ftr: luminance8 or angle4 or luminance3; noi: noise
        for nz = 1:2
          Fm{s,m,i,nz} = [hiss{s,m}(i+3*(nz-1)).his.ftr(1,:)' hiss{s,m}(i+3*(nz-1)).his.qu1' ...
              hiss{s,m}(i+3*(nz-1)).his.qu2' qu12' (nz-1)*ones(nt1, 1)];
          Ft_smi{s,m,i,nz} = array2table(Fm{s,m,i,nz}, 'VariableNames', dsf);
        end
        %Fm{s,m,i} = [hiss{s,m}(i).his.ftr(1,:)' hiss{s,m}(i).his.qu1' ...
        %    hiss{s,m}(i).his.qu2' qu12' zeros(nt0, 1)]; 
        %Fm1{s,m,i} = [hiss{s,m}(i+3).his.ftr(1,:)' hiss{s,m}(i+3).his.qu1' ...
        %    hiss{s,m}(i+3).his.qu2' qu12n1' ones(nt1, 1)];
        %Ft_smi1{s,m,i} = array2table(Fm1{s,m,i}, 'VariableNames', dsf);
        %Ft_smi{s,m,i} = array2table(Fm{s,m,i}, 'VariableNames', dsf);
      end      
      
    end % i
  end % m
end % s 