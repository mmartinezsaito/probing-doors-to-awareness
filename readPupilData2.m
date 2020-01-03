% convert P M i_events
%  for E3_4: into cells with 3 elements (detec, discr1, discr2)

datapath2 = '~/Data/PupilSize/';
datapath4 = '/home/mario/MEGA/Neuroscience/Data/Pupillometry/E3_4/';	
datapath7 = '/home/mario/MEGA/Neuroscience/Data/Pupillometry/E3_7/';	

for s = 7%1:10
  for m = 3%  [1 3]
    fprintf(1, '\n\nsubject %d, measure %d\n', s, m)
    if s <= 6 % E3.4
      load([datapath4 's' sprintf('%02d', s) 'm' sprintf('%1d', m) '.mat']) 

      if s == 1 && m == 1 % ralfCR
        ie23 = find(ismember(i_events{2},[1 2])); 
        i_events{2}([1:ie23(4)]) = []; i_events{3}([1:ie23(4)]) = []; 
        it1 = [P{2}{:,1}] > M{2}{ie23(4), 1};
        M{2}(1:ie23(4), :) = []; 
        P{2} = P{2}(it1, :); 
      elseif s == 6 && m == 3  % saadatPAS
        ie1 = find(ismember(i_events{1},[1 2])); 
        i_events{1} = i_events{1}(ie1(2):end); 
      end
      % (1) detec: 720; (2) discr1: 360; (3) discr2: 270
      if s == 2 && m == 3 % kikaPAS
        ie1 = find(ismember(i_events{1},[1 2])); ie23 = find(ismember(i_events{2},[1 2])); 
        i_events{2} = i_events{1}(796:end); %ie1(5)=796
        i_events{1}(796:end) = [];
        i_events{2} = [i_events{2}; i_events{3}(1:241)];  % ie23(1)=241
        i_events{3}(1:241)= [];
        it12 = [P{1}{:,1}] < M{1}{796, 1};
        it23 = [P{2}{:,1}] < M{2}{241, 1};
        M{3} = M{2};                   P{3} = P{2};
        M{2} = M{1}(796:end, :);       
        M{1}(796:end,:) = [];          
        M{2} = [M{2}; M{3}(1:241, :)];
        M{3}(1:241,:)= [];             
        P{2} = P{1}(it12, :); 
        P{1}(it12,:) = [];
        P{2} = [P{2}; P{3}(it23, :)];
        P{3}(it23,:)= [];
      else
        ie23 = find(ismember(i_events{2}, [1 2])); 
        i_events{2}(ie23(5):end) = []; 
        i_events{3}(1:ie23(4)) = []; 
        it2 = [P{2}{:,1}] < M{2}{ie23(5), 1};
        it3 = [P{2}{:,1}] > M{2}{ie23(4), 1};
        M{3} = M{2}(ie23(5):end, :); 
        M{2} = M{2}(1:ie23(4), :); 
        P{3} = P{2}(it3, :);
        P{2} = P{2}(it2, :);
      end

      for i = 2; % 1:3
      end

    elseif s > 6 % E3.7
      load([datapath7 's' sprintf('%02d', s) 'm' sprintf('%1d', m) '.mat']) 

      if s == 7 && m == 3  % feyzaPASn1
        for i = 4:6, i_events{i} = [1; i_events{i}]; end
        for i = 2:6 , M{i}{1, 2} = 't0'; end
      elseif s == 8 && m == 1  % barriCRn0
        i_events{1}(336:end) = []; i_events{2}(1:32) = []; 
        i_events{1} = [i_events{1}; i_events{2}]; for i=2:3, i_events{i} = i_events{1}; end
        it01 = [P{1}{:,1}] >= M{1}{336, 1}; it02 = [P{2}{:,1}] <= M{2}{32, 1};
        M{1}(336:end,:) = []; M{2}(1:32,:) = []; M{1} = [M{1}; M{2}]; M{2} = M{3}; 
        P{1}(it01,:) = []; P{2}(it02,:) = []; P{1} = [P{1}; P{2}]; P{2} = P{3};
      elseif s == 8 && m == 3   % barriCRn1
        for i = 1:3, i_events{i} = [1; i_events{i}]; M{i}{1, 2} = 't0'; end
      elseif s ==10 && m == 1   % himajinCRn0
        for i = 4:6, i_events{i} = [1; i_events{i}]; M{i}{1, 2} = 't0'; end
      end
      if 1
        ie0 = find(ismember(i_events{1}, [1 2])); ie1 = find(ismember(i_events{4}, [1 2])); 
        i_events{1}(ie0(5):end) = [];             i_events{4}(ie1(5):end) = []; 
        i_events{2}([1:ie0(4) ie0(9):end]) = [];  i_events{5}([1:ie1(4) ie1(9):end]) = []; 
        i_events{3}(1:ie0(8)) = [];               i_events{6}(1:ie1(8)) = []; 
        it01 = [P{1}{:,1}] <= M{1}{ie0(4), 1};                                  it11 = [P{2}{:,1}] <= M{2}{ie1(4), 1};
        it02 = [P{1}{:,1}] >= M{1}{ie0(5), 1} & [P{1}{:,1}] <= M{1}{ie0(8), 1}; it12 = [P{2}{:,1}] >= M{2}{ie1(5), 1} & [P{2}{:,1}] <= M{2}{ie1(8), 1};
        it03 = [P{1}{:,1}] >= M{1}{ie0(9), 1};                                  it13 = [P{2}{:,1}] >= M{2}{ie1(9), 1};
        M{4} = M{2};
        M{3} = M{1}(ie0(9):end, :);    M{6} = M{4}(ie1(9):end, :); 
        M{2} = M{1}(ie0(5):ie0(8), :); M{5} = M{4}(ie1(5):ie1(8), :); 
        M{1}(ie0(5):end, :) = [];      M{4}(ie1(5):end, :) = []; 
        P{4} = P{2};
        P{3} = P{1}(it03, :); P{6} = P{4}(it13, :); 
        P{2} = P{1}(it02, :); P{5} = P{4}(it12, :);
        P{1} = P{1}(it01, :); P{4} = P{4}(it11, :);
      end
      i = 1; % 1:6
    end % E3.4 or E3.7
    save(['~/Desktop/s' sprintf('%02d', s) 'm' num2str(m) 'c.mat'], 'P', 'M', 'i_events') 
  end % m  
end % s

%{
if s <= 6 
elseif s > 6
  if tfpc(1) == 3 % for feyzaPASn1
    i_t0 = [0; i_t0];
  elseif all(tfpc == [4 2]') || all(tfpc == [4 3]') % for barriCRn0
    i_t0 = [0; i_t0];
  elseif tfpc(1) == 5  % for barriCRn1
    i_t0 = [0; i_t0];
  elseif tfpc(1) == 6  % for himajinCRn0
    i_t0 = [0; i_t0];
  end 
  if all(tfpc == [4 2]') || all(tfpc == [4 3]') % for barriCRn0 
    ntchk = i_tend(2*(tfpc(2)-1)) - i_t0(2*(tfpc(2)-1)) - 1;  
  else      
    ntchk = i_tend(2*tfpc(2)) - i_t0(2*tfpc(2)) - 1;  
  end
  if all(tfpc == [4 1]') % for barriCRn0
    %fprintf(1, '\nfaulty file: 2nd t0 starts at %d', i_t0(4)) 
  end  
end
%}
