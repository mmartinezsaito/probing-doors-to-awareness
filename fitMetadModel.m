 % fill fit{} with needed subjects data
  nR_S1_adj = {}; nR_S2_adj = {}; fit = {}; fit_adj = {};
  for s = ss
    if strcmp(subgrp, 'avg') || strcmp(subgrp, 'all'), si = s; end
    
    for l = lvs
      if islr && tk==2
        lr = Ft{any(Ft.sid==si, 2) & Ft.task==tk & Ft.mjdg==mj, 'ftr1_raw'} > 4;
        q2_l = Ft{any(Ft.sid==si, 2) & Ft.task==tk & Ft.ftr1_raw==5-l & Ft.mjdg==mj, 'qu2'};
        q2_r = Ft{any(Ft.sid==si, 2) & Ft.task==tk & Ft.ftr1_raw==4+l & Ft.mjdg==mj, 'qu2'}; 
        q1_l = ~Ft{any(Ft.sid==si, 2) & Ft.task==tk & Ft.ftr1_raw==5-l & Ft.mjdg==mj, 'qu1'};
        q1_r = Ft{any(Ft.sid==si, 2) & Ft.task==tk & Ft.ftr1_raw==4+l & Ft.mjdg==mj, 'qu1'}; 
        nR_S1 = sum((2 * q1_l - 1) .* q2_l == uq12);     
        nR_S2 = sum((2 * q1_r - 1) .* q2_r == uq12);     
      elseif tk==1 && mj==3 && ist1uv
        q2S1 = Ft{any(Ft.sid==si, 2) & Ft.task==tk & Ft.ftr1==1 & Ft.mjdg==mj, 'qu2'}; 
        q2S2 = Ft{any(Ft.sid==si, 2) & Ft.task==tk & Ft.ftr1==l & Ft.mjdg==mj, 'qu2'}; 
        % 2 points  
        %nR_S1 = sum((q2S1~=1) == uq1);     
        %nR_S2 = sum((q2S2~=1) == uq1);     
        % 4 points  
        nR_S1 = sum(q2S1 == uq2);     
        nR_S2 = sum(q2S2 == uq2);    
      else
        nR_S1 = sum(Ft{any(Ft.sid==si, 2) & Ft.task==tk & Ft.ftr1==1 & Ft.mjdg==mj, 'qu12'} == uq12); 
        nR_S2 = sum(Ft{any(Ft.sid==si, 2) & Ft.task==tk & Ft.ftr1==l & Ft.mjdg==mj, 'qu12'} == uq12);   
        
        if tk==1 && mj==1 && ist1uv && iscr2vr  % collapse -4:-1 in CR to 1 in VR
          nR_S1 = cr2vr(nR_S1); nR_S2 = cr2vr(nR_S2);
        end
      end
      
      if ist2uv && (mj~=3 || tk~=1) && is4
        s1divs2 = G4f13{G4f13.mjdg==mj_s&G4f13.task==tk&G4f13.sid==s&G4f13.ftr1==l, 's'};
      end
     
      % ugly patch
      if mj==1 && tk==1 && ~ist1uv
        if     s==3 && ~strcmp(subgrp, 'pool') 
          switch l
          case 2,  nR_S2 = nR_S2 + [-1  0  0 0 1 0 0 0];
          case 3,  nR_S2 = nR_S2 + [-2 -1  0 0 2 1 0 0];
          case 4,  nR_S2 = nR_S2 + [-2 -1 -1 0 2 1 1 0]; end
        elseif s==4 && ~strcmp(subgrp, 'pool') 
          if l==2, nR_S2 = nR_S2 + [-1  0  0 0 1 0 0 0]; end
        elseif s==5 && ~strcmp(subgrp, 'pool')
          if l==4, nR_S2 = nR_S2 + [0  0 -10 0 5 3 2 0]; end
        end
      end
    
      adj_f = 1/length(nR_S1);
      nR_S1_adj{s} = nR_S1 + adj_f; nR_S2_adj{s} = nR_S2 + adj_f;
      try
        if isempty(rstr) && ~ist1uv, fit_adj{s}(l) = fit_meta_d_MLE(nR_S1_adj{s}, nR_S2_adj{s}, s1divs2);    
        elseif ~ist1uv,              fit_adj{s}(l) = fit_rs_meta_d_MLE(nR_S1_adj{s}, nR_S2_adj{s}, s1divs2);          
        elseif ist1uv,               fit_adj{s}(l) = type1_SDT(nR_S1_adj{s}, nR_S2_adj{s}, iseqvar, nR_thr);          
        end
      catch ME
        if strcmp(ME.identifier, 'optim:barrier:UsrObjUndefAtX0')
          warning('adding noise to nRs')
          if isempty(rstr) && ~ist1uv, fit_adj{s}(l) = fit_meta_d_MLE(nR_S1_adj{s}+adj_f*rand(1,length(nR_S1)), nR_S2_adj{s}+adj_f*rand(1,length(nR_S2)), s1divs2); 
          elseif ~ist1uv,              fit_adj{s}(l) = fit_rs_meta_d_MLE(nR_S1_adj{s}+adj_f*rand(1,length(nR_S1)), nR_S2_adj{s}+adj_f*rand(1,length(nR_S2)), s1divs2);
          elseif ist1uv,               fit_adj{s}(l) = type1_SDT(nR_S1_adj{s}+adj_f*rand(1,length(nR_S1)), nR_S2_adj{s}+adj_f*rand(1,length(nR_S2)), iseqvar, nR_thr);          
          end
        else, rethrow(ME)
        end       
      end
      fit{s} = fit_adj{s};
      
      %G{G.sid==s & G.task==tk & G.mjdg==mj & G.ftr1==l & Gf2ni,'da'} = fit{s}(l).da;
      G{G.sid==s & G.task==tk & G.mjdg==mj & G.ftr1==l,'da'} = fit{s}(l).da;
      if ~ist1uv
        G{G.sid==s & G.task==tk & G.mjdg==mj & G.ftr1==l,['meta_da' rstr]} = getfield(fit{s}(l), ['meta_da' rstr]); 
        if isempty(rstr), t1ORmeta = 'meta_ca'; else, t1ORmeta = 't1ca'; end 
        G{G.sid==s & G.task==tk & G.mjdg==mj & G.ftr1==l,['t1ca' rstr]} = getfield(fit{s}(l), [t1ORmeta rstr]);  
      else
        G{G.sid==s & G.task==tk & G.mjdg==mj & G.ftr1==l, 's'} = fit{s}(l).s;
        G{G.sid==s & G.task==tk & G.mjdg==mj & G.ftr1==l, 'c_1'} = fit{s}(l).c_1;
        G{G.sid==s & G.task==tk & G.mjdg==mj & G.ftr1==l, 'd_1'} = fit{s}(l).d_1;
      end
    end  
    
    if strcmp(subgrp, 'one') || strcmp(subgrp, 'pool')
      break
      % fit{1} for one, fit{1} for pool
    elseif s >= ss(end) && strcmp(subgrp, 'avg') 
      fnc = fieldnames(fit{ss(1)});
      for l = lvs
        fit{ss(end)+1}(l) = fit{ss(1)}(1); fit{ss(end)+2}(l) = fit{ss(1)}(1);
        fni = 1:8;
        if ~isempty(rstr), fni = [fni 14:18]; 
        elseif ist1uv,     fni = 1:5;
        end
        for i = fni %length(fnc)
          tmp = {}; for s = ss, tmp{s} = fit{s}(l).(fnc{i}); end
          if (~ist1uv && isempty(rstr) && ismember(i,[7 8])) || ...
             (~ist1uv && ~isempty(rstr) && ismember(i,[6 16])) || ...
             ( ist1uv && i==5), tmp = cell2mat(tmp');
          else,                 tmp = [tmp{:}];  
          end
          
          [mn, se] = grpstats(tmp, [], {'mean', 'sem'});
          fit{ss(end)+1}(l).(fnc{i}) = mn; fit{ss(end)+2}(l).(fnc{i}) = se; 
        end
      end
    end    
  end