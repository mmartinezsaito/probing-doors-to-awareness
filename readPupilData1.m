function [P, M, i_events, pars, hiss] = readPupilData1(varargin)

datapath = '/home/mario/MEGA/Neuroscience/Data/Pupillometry/';	

loadDataFilenames

txtfpartsE3_4 = [1 2 2; 1 1 2];
txtfpartsE3_7 = [1 1 1 2 2 2; 1 2 3 1 2 3];

mregex = 'E3_\d(?<dedi>\w{5}\d?)(?<meta>[A-Z]{2,3})_(?<date>\d{8})T(?<hour>\d{6})\.mat';
tregex = 'E3_\dD(?<date>\d{6})(?<sid>[a-z]+)(?<meta>[A-Z]{2,3})_?(?<part>.*?)\.txt';

if nargin == 5 % txtf1 txtf2 matf1 matf2 matf3
  txtf1 = varargin{1};
  txtf2 = varargin{2};
  %Sp = regexp(txtf1, tregex, 'tokens');
  %Sp_split = regexp(txtf1, '[12]of2', 'split');
  P = {}; M = {};
  [P{1}, M{1}] = readtxtf(txtf1);
  [P{2}, M{2}] = readtxtf(txtf2);
  matf = {}; i_events = {};
  txtfparts = [1 2 2; 1 1 2];
  for i = 1:3
    matf{i} = varargin{i + 2};
    pars(i) = load(matf{i}, 'par'); 
    hiss(i) = load(matf{i}, 'his'); 
    i_events{i} = chknt(M{txtfparts(1, i)}, txtfparts(:, i), pars(i).par);
  end
elseif nargin == 0
  %matfl = dir('*.mat'); 
  %Sm = regexp([matfl.name], mregex, 'tokens'); %'names');
  %txtfl = dir('*.txt'); 
  %Sp = regexp([txtfl.name], tregex, 'tokens');
  %dateSm = cellfun(@(x) x(3), Sm);
  %dateSp = cellfun(@(x) x(1), Sp);
  for s = 2%1:10
    switch s > 6, case 0, cd([datapath 'E3_4']), case 1, cd([datapath 'E3_7']), end
    for m = 3% [1 3]
      fprintf(1, '\n\nsubject %d, measure %d\n', s, m)
      P = {}; M = {}; i_events = {}; i_t0 = {}; i_onsetst = {}; i_tend = {};
      if s <= 6 % E3.4
        tfpm = txtfpartsE3_4(1,:);
        tfp = txtfpartsE3_4;
        if s == 1 && m == 1, tfp = [1 3 3; 1 1 2]; end  % ralfCR
        if s == 2 && m == 3, tfpm = [1 1 2]; tfp = [1 4 4; 1 1 2]; end  % kikaPAS
        if s == 6 && m == 3, tfp = [5 2 2; 1 1 2]; end  % saadatPAS
        for i = 1:2 
          [P{i}, M{i}] = readtxtf(sfn{s, m}{i}); 
        end 
        for i = 1:3
          pars(i) = load(sfn{s, m + 1}{i}, 'par');
          hiss(i) = load(sfn{s, m + 1}{i}, 'his');
          [i_events{i}, i_t0{i}, i_onsetst{i}, i_tend{i}] = chknt(M{tfpm(i)}, tfp(:, i), pars(i).par);
        end
      elseif s >= 7 % E3.7
        tfpm = txtfpartsE3_7(1, :); 
        tfp = txtfpartsE3_7;
        if s == 7 && m == 3, tfp(1, :) = [1 1 1 3 3 3]; end  % feyzaPASn1
        if s == 8 && m == 1, tfpm = [1 2 2 3 3 3]; tfp(1, :) = [4 4 4 2 2 2]; [P{3}, M{3}] = readtxtf(sfn{s, m}{3}); end  % barriCRn0
        if s == 8 && m == 3, tfp(1, :) = [5 5 5 2 2 2]; end  % barriCRn1
        if s ==10 && m == 1, tfp(1, :) = [1 1 1 6 6 6]; end  % himajinCRn0
        for i = 1:2 
          [P{i}, M{i}] = readtxtf(sfn{s, m}{i}); 
        end 
        for i = 1:6
          pars(i) = load(sfn{s, m + 1}{i}, 'par');
          hiss(i) = load(sfn{s, m + 1}{i}, 'his');
          [i_events{i}, i_t0{i}, i_onsetst{i}, i_tend{i}] = chknt(M{tfpm(i)}, tfp(:, i), pars(i).par);
        end
      end % E3.4 or E3.7
      save(['~/Desktop/s' sprintf('%02d', s) 'm' num2str(m) '.mat']) 
    end % m  
  end % s
end



function [P, M] = readtxtf(txtf)
% extract columns for P: Time[ms] LDia[cpx] LArea{c(px2)] LPORX[cpx] LPORY[cpx]
% extract columns for M: Time[ms] Event 

fs = 240; decf = 1; 

fid = fopen(txtf, 'r');

for i = 1:100
  tline = fgets(fid); 
  if strncmp(tline, ['Time' char(9) 'Type'], 9)
    break
  end
  if i == 100, error('Wrong file?'); end
end 

P = {}; M = {};
p = 0; m = 0; s = 0;
while 1
  tline = fgetl(fid);
  if tline == -1
    break
  elseif strfind(tline, '# Message:');
    m = m + 1;
    l = sscanf(tline, '%f%*s%*u%*s%*s%s'); 
    M{m, 1} = uint32(l(1) / 10^3); % convert micros to milis 
    M{m, 2} = char(l(2:end))'; 
  else
    s = s + 1;
    if mod(s, decf) == 0 % decimation factor 10: downsample from 240Hz to 24Hz 
      % textscan(fid, '%u%*s%u%f%f%f%f%*d', 'CollectOutput', 1, 'CommentStyle', '#', 'HeaderLines', 22);
      p = p + 1;
      l = sscanf(tline, '%f%*s%*u%f%f%f%f%*d');
      P{p, 1} = uint32(l(1) / 10^3); % convert micros to milis 
      P{p, 2} = uint16(l(2) * 100); % convert to cpx 
      P{p, 3} = uint32(l(3) * 100); % convert to c(px^2)
      P{p, 4} = uint32(l(4) * 100); % convert to cpx
      P{p, 5} = uint32(l(5) * 100); % convert to cpx
    end
  end
end
fseek(fid, 0, 'bof');
fclose(fid);
fprintf(1, '%d messages, %d extracted data lines, %d read data lines\n', m, p, s)



function [i_events, i_t0, i_onsetst, i_tend] = chknt(M, tfpc, par)

A1 = cellfun(@(x) strcmp(x, 't0'),  M(:,2));
i_t0 = find(A1);
A2 = cellfun(@(x) strcmp(x, 't_end'),  M(:,2));
i_tend = find(A2);
A0 = cellfun(@(x) strcmp(x, 'onset_st'),  M(:,2));
i_onsetst = find(A0);
i_events = A1 + 2 * A2;

if strcmp(par.filename(1:4), 'E3_4')
  if tfpc(1) == 1
    ntchk = i_tend(2) - i_t0(2) - 1;  
  elseif tfpc(1) == 2
    ntchk = i_tend(tfpc(2)*2) - i_t0(tfpc(2)*2) - 1;  
  elseif tfpc(1) == 3 % for ralfCR
    ntchk = i_tend(tfpc(2)*2+2) - i_t0(tfpc(2)*2+2) - 1;  
  elseif tfpc(1) == 4 % for kikaPAS discr
    if tfpc(2) == 1
      ntchk = size(M,1) - i_t0(4) - 1; % 950-837-1=112 
    elseif tfpc(2) == 2   
      ntchk = i_tend(3) - i_t0(2) - 1;  
      ntchk2 = i_tend(1); % 241 (from discr1) 
    end %241+112=353
  elseif tfpc(1) == 5 % for saadatPAS
    ntchk = i_tend(2) - i_t0(3) - 1;  
  end
  if all(tfpc == [4 1]')
    fprintf(1, '\nfaulty file: %d trials', ntchk) 
  elseif all(tfpc == [4 2]')
    fprintf(1, '\n additional %d trials from previous faulty file', ntchk2) 
  end
elseif strcmp(par.filename(1:4), 'E3_7')
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

if par.nt ~= ntchk
  warning('Number of trials don''t match: %d ~= %d', par.nt, ntchk)
else
  fprintf(1, '\nnt = %d', par.nt)
end


