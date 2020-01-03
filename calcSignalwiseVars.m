 q1={}; q2={}; q12={}; idx1={}; jdx1={}; idx2={}; jdx2={}; idx12={}; jdx12={};
  labels1={}; labels2={}; labels12={}; labels32={}; scores1={}; scores2={}; scores12={};  
  nwq1 = linspace(0, 1, 2); nwq12 = linspace(0, 1, 8); nwq2 = linspace(0, 1, 4);
  for l = lvs
    lidx = any(Ft.ftr1==[1 l], 2);  % lvs
       
    for s = ss
      if strcmp(subgrp, 'avg') || strcmp(subgrp, 'all'), si = s; end  
      Ft_stml = Ft(any(Ft.sid==si, 2) & Ft.task==tk & Ft.mjdg==mj & lidx, :)
      
      ft1 = Ft_stml.ftr1;
      % perceptual decisions
      if islr && tk==2
        lr = Ft_stml.ftr1_raw > 4;
        q1{s,l} = ~xor(lr, Ft_stml.qu1); q2{s,l} = Ft_stml.qu2; 
        q12{s,l} = (2 * q1{s,l} - 1) .* q2{s,l};  
      else
        q1{s,l} = Ft_stml.qu1; q2{s,l} = Ft_stml.qu2; q12{s,l} = Ft_stml.qu12;  
      end
            
      % ROC curves with Pointwise Confidence Bounds      
      [idx1{s,l} jdx1{s,l}] = find(q1{s,l}' == uq1');
      [idx2{s,l} jdx2{s,l}] = find(q2{s,l}' == uq2');
      [idx12{s,l} jdx12{s,l}] = find(q12{s,l}' == uq12');
      if mj==3 && tk==1   % reconstruct a new qu1 for PAS
        q1{s,l} = q2{s,l}~=1; idx1{s,l} = (idx2{s,l}~=1)+1;
      end
    
      % Type 1 condition: nothing(1) versus something(2:end)
      if islr && tk==2, labels1{s,l} = Ft_stml.ftr1_raw > 4;
      else,             labels1{s,l} = ft1(jdx1{s,l})~=1;
      end
      if any(isnan(q2{s,l})), ih = find(isnan(q2{s,l})); q1{s,l}(ih) = []; end    
      scores1{s,l} = nwq1(idx1{s,l})';
            
      % Type 2 condition: correct qu1 answers for level l of stimulus
      if tk==1
        labels2{s,l} = (labels1{s,l} & q1{s,l}) | (~labels1{s,l} & ~q1{s,l}); % accuracy
      else
        labels2{s,l} = Ft_stml.qu1;       
      end
      scores12{s,l} = nwq12(idx12{s,l})';
      if any(isnan(q2{s,l}))
        ih = find(isnan(q2{s,l}));
        [~, ind] = sort([setdiff(1:length(q2{s,l}), ih) ih']);
        tmp = [scores12{s,l}' NaN(size(scores12{s,l},2), length(ih))]'; 
        scores12{s,l} = tmp(ind);
        %scores12{s,l} = [scores12{s,l}(1:ih-1); NaN; scores12{s,l}(ih:end)]; 
      end
      scores2{s,l} = nwq2(idx2{s,l})';

      %labels3{s,l} = (labels2{s,l} & q2{s,l}>2) | (~labels2{s,l} & q2{s,l}<3); 
      labels32{s,l} = ((2*labels2{s,l}-1) .* (q2{s,l}-2.5) + 1.5)/3; % accuracy
      
      if strcmp(subgrp, 'pool') || strcmp(subgrp, 'one'), break, end      
    end    
  end
  Ft_tm = Ft(Ft.task==tk & Ft.mjdg==mj, :);
  if islr && tk==2
    lr = Ft_tm.ftr1_raw > 4;
    Ft_tm.qu1 = ~xor(lr, Ft_tm.qu1); Ft_tm.qu12 = (2 * Ft_tm.qu1 - 1) .* Ft_tm.qu2;  
  end
  
  % subject-wise averaged
  q1M = reshape(grpstats(Ft_tm.qu1, {Ft_tm.ftr1 Ft_tm.sid}), length(ss), []);
  q1Mrm = grpstats(Ft_tm.right1, {Ft_tm.ftr1 Ft_tm.sid}); ct = crosstab(Ft_tm.ftr1,Ft_tm.sid)'; ct(find(ct)) = q1Mrm(:); ct(find(ct==0)) = NaN; q1Mr = ct; % accuracy 
  q1Mhm = grpstats(Ft_tm.hit, {Ft_tm.ftr1 Ft_tm.sid}); ct = crosstab(Ft_tm.ftr1,Ft_tm.sid)'; ct(find(ct)) = q1Mhm(:); ct(find(ct==0)) = NaN; q1Mh = ct; % hits 
  q1Mcm = grpstats(Ft_tm.correj, {Ft_tm.ftr1 Ft_tm.sid}); ct = crosstab(Ft_tm.ftr1,Ft_tm.sid)'; ct(find(ct)) = q1Mcm(:); ct(find(ct==0)) = NaN; q1Mc = ct; % correj 
  q2M = reshape(grpstats((Ft_tm.qu2-1)/3, {Ft_tm.ftr1 Ft_tm.sid}), length(ss), []);
  q2M0m = grpstats((Ft_tm{Ft_tm.qu1==0,'qu2'}-1)/3, {Ft_tm{Ft_tm.qu1==0,'ftr1'} Ft_tm{Ft_tm.qu1==0,'sid'}}); ct = crosstab(Ft_tm{Ft_tm.qu1==0,'ftr1'},Ft_tm{Ft_tm.qu1==0,'sid'})'; ct(find(ct)) = q2M0m(:); ct(find(ct==0)) = NaN; q2M0 = ct; % no
  q2M1m = grpstats((Ft_tm{Ft_tm.qu1==1,'qu2'}-1)/3, {Ft_tm{Ft_tm.qu1==1,'ftr1'} Ft_tm{Ft_tm.qu1==1,'sid'}}); ct = crosstab(Ft_tm{Ft_tm.qu1==1,'ftr1'},Ft_tm{Ft_tm.qu1==1,'sid'})'; ct(find(ct)) = q2M1m(:); ct(find(ct==0)) = NaN; q2M1 = ct; % yes  
  q2Mwm = grpstats((Ft_tm{Ft_tm.wrong1,'qu2'}-1)/3, {Ft_tm{Ft_tm.wrong1,'ftr1'} Ft_tm{Ft_tm.wrong1,'sid'}}); ct = crosstab(Ft_tm{Ft_tm.wrong1,'ftr1'},Ft_tm{Ft_tm.wrong1,'sid'})'; ct(find(ct)) = q2Mwm(:); ct(find(ct==0)) = NaN; q2Mw = ct; % wrong1
  q2Mrm = grpstats((Ft_tm{Ft_tm.right1,'qu2'}-1)/3, {Ft_tm{Ft_tm.right1,'ftr1'} Ft_tm{Ft_tm.right1,'sid'}}); ct = crosstab(Ft_tm{Ft_tm.right1,'ftr1'},Ft_tm{Ft_tm.right1,'sid'})'; ct(find(ct)) = q2Mrm(:); ct(find(ct==0)) = NaN; q2Mr = ct; % right1 
  q21M = reshape(grpstats(Ft_tm.qu2>1, {Ft_tm.ftr1 Ft_tm.sid}), length(ss), []);
  for i=1:4, q2Mrating{i}=reshape(grpstats(Ft_tm.qu2==i, {Ft_tm.ftr1 Ft_tm.sid}), length(ss), []); end
  %for i=1:4, q2M0rating{i}=reshape(grpstats(Ft_tm{Ft_tm.qu1==0,'qu2'}==i, {Ft_tm{Ft_tm.qu1==0,'ftr1'} Ft_tm{Ft_tm.qu1==0,'sid'}}), length(ss), []); end
  %for i=1:4, q2M1rating{i}=reshape(grpstats(Ft_tm{Ft_tm.qu1==1,'qu2'}==i, {Ft_tm{Ft_tm.qu1==1,'ftr1'} Ft_tm{Ft_tm.qu1==1,'sid'}}), length(ss), []); end
  for i=1:4, q2M0rating{i}=grpstats(Ft_tm{Ft_tm.qu1==0,'qu2'}==i, {Ft_tm{Ft_tm.qu1==0, 'sid'} Ft_tm{Ft_tm.qu1==0,'ftr1'}}); ct=crosstab(Ft_tm{Ft_tm.qu1==0,'ftr1'},Ft_tm{Ft_tm.qu1==0,'sid'})'; ct(find(ct))=q2M0rating{i}(:); ct(find(ct==0))=NaN; q2M0rating{i}=ct; end
  for i=1:4, q2M1rating{i}=grpstats(Ft_tm{Ft_tm.qu1==1,'qu2'}==i, {Ft_tm{Ft_tm.qu1==1, 'sid'} Ft_tm{Ft_tm.qu1==1,'ftr1'}}); ct=crosstab(Ft_tm{Ft_tm.qu1==1,'ftr1'},Ft_tm{Ft_tm.qu1==1,'sid'})'; ct(find(ct))=q2M1rating{i}(:); ct(find(ct==0))=NaN; q2M1rating{i}=ct; end
  % pooled 
  [q1mp, q1sp] = grpstats(Ft_tm.qu1, {Ft_tm.ftr1}, {'mean', 'sem'});
  [q2mp, q2sp] = grpstats((Ft_tm.qu2-1)/3, {Ft_tm.ftr1}, {'mean', 'sem'}); 
  [q2mp0, q2sp0] = grpstats((Ft_tm{Ft_tm.qu1==0,'qu2'}-1)/3, {Ft_tm{Ft_tm.qu1==0,'ftr1'}}, {'mean', 'sem'}); 
  [q2mp1, q2sp1] = grpstats((Ft_tm{Ft_tm.qu1==1,'qu2'}-1)/3, {Ft_tm{Ft_tm.qu1==1,'ftr1'}}, {'mean', 'sem'}); 
  [q21mp, q21sp] = grpstats(Ft_tm.qu2>1, {Ft_tm.ftr1}, {'mean', 'sem'});