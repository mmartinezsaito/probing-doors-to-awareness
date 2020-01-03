% ATTENTION TO PATCH
% Ft7{Ft7.task==2,'qu1'} = double(~Ft7{Ft7.task==2,'qu1'})
% Ft7{Ft7.task==2,'qu12'} = -Ft7{Ft7.task==2,'qu12'}

Ft = {}; 
if length(ss) == 6 && all(ss == 1:6)
  for m = ms
    for i = is
      Ft = [Ft; table(m*ones(size(Ft_mi{m,i},1),1), i*ones(size(Ft_mi{m,i},1),1), ...
          'VariableNames', {'mjdg', 'task'}) Ft_mi{m, i}]; 
    end
  end, Ft.Properties.VariableNames

  % renaming ftr1 of discr1
  olb = Ft{Ft.task==2, 'ftr1'}; nlb = olb;
  for i = 1:4, nlb(any(olb==[i 9-i], 2)) = 5-i; end
  Ft{Ft.task==2,'ftr'} = olb; 
  Ft{Ft.task==2,'ftr1'} = nlb;
   
  % making ftr2 categorical
  Ft.ftr2 = categorical(Ft.ftr2, 'Ordinal', false);
else
  for m = ms
    for i = is
      Ft = [Ft; table(m*ones(size(Ft_mi{m,i},1),1), i*ones(size(Ft_mi{m,i},1),1), ...
        'VariableNames', {'mjdg', 'task'}) Ft_mi{m, i}]; 
    end
  end, Ft.Properties.VariableNames
  
  %renaming ftr to ftr1
  Ft.Properties.VariableNames{6} = 'ftr1';
  
  % making noi categorical
  Ft.noi = categorical(Ft.noi, 'Ordinal', false);
end

% extracting pupil size average by clipping with a20s window
cnv = size(Ft, 2);
for ei = 1:7
  Ft{:, cnv+ei} = zeros(size(Ft,1), 1); % make new column
  Ft.Properties.VariableNames{cnv+ei} = ['ps_' lkev{ei}]; % name it  
  for s = ss
    for m = ms
      for i = is
          
        if s < 7    
          A = mean(W_smie{s,m,i,ei}(:, a20s), 2); 
          try
            Ft{Ft.task==i & Ft.mjdg==m & Ft.sid==s, cnv+ei} = A;    
          catch ME
            A = [A; NaN(sum(Ft.task==i & Ft.mjdg==m & Ft.sid==s)-length(A), 1)];    
            Ft{Ft.task==i & Ft.mjdg==m & Ft.sid==s, cnv+ei} = A;
          end
        else
          for nz = 1:2  
            A = mean(W_smie{s,m,i,ei,nz}(:, a20s), 2); 
            try
              Ft{Ft.task==i & Ft.mjdg==m & Ft.sid==s & Ft.noi==categorical(nz-1), cnv+ei} = A;    
            catch ME
              A = [A; NaN(sum(Ft.task==i & Ft.mjdg==m & Ft.sid==s & Ft.noi==categorical(nz-1))-length(A), 1)];    
              Ft{Ft.task==i & Ft.mjdg==m & Ft.sid==s & Ft.noi==categorical(nz-1), cnv+ei} = A;
            end
          end % nz                
        end 
        
      end % i
    end % m
  end % s 
end % ei   

% removing PS outliers (beyond 3-sigma)
ps3sigma = std(Ft{~any(isnan(Ft{:, cnv+(1:ei)}), 2), cnv+(1:ei)}, 1) * 3; 
ps3sigma_avg = mean(ps3sigma);
outlier_rows = any(Ft{:,cnv+(1:ei)} > ps3sigma_avg, 2); % 3.4: 543 / 16200 = 3.35%
%Ft = Ft(~outlier_rows, :);

% removing missing values 
missing_rows = any(ismissing(Ft), 2);
Ft = Ft(~missing_rows, :);

% computing confusion matrix variables
Ft.hit = Ft.ftr1~=1 & Ft.qu1;    
Ft.falsal = Ft.ftr1==1 & Ft.qu1;
Ft.miss = Ft.ftr1~=1 & ~Ft.qu1;
Ft.correj = Ft.ftr1==1 & ~Ft.qu1;
% ATTENTION: THESE 4 ABOVE ARE MEANINGLESS FOR THE TASK 2
Ft.right1 = Ft.hit | Ft.correj;
Ft.wrong1 = Ft.miss | Ft.falsal;
[~, jdx] = find(Ft{:, {'hit' 'falsal' 'miss' 'correj'}}); Ft.cm = jdx;

Ft.correctsure = Ft.right1 & Ft.qu2>2.5;
Ft.overconf = Ft.wrong1 & Ft.qu2>2.5;
Ft.underconf = Ft.right1 & Ft.qu2<2.5;
Ft.correcthesitant = Ft.wrong1 & Ft.qu2<2.5;
Ft.right2 = Ft.correctsure | Ft.correcthesitant;
Ft.wrong2 = Ft.underconf | Ft.overconf;