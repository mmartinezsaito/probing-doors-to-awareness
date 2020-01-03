if topool  
  Ft_mi = {}; Fma = {}; rti_qu1a = {}; rti_qu2a = {};    
  for i = is
    for m = ms
      Ft_mi{m,i}=[]; Fma{m,i}=[]; 
      rti_qu1a{m,i} = []; rti_qu2a{m,i} = [];
      for s = ss
        if    s <= 6
          rti_qu1a{m,i} = [rti_qu1a{m,i} rti_qu1{s,m,i}]; 
          rti_qu2a{m,i} = [rti_qu2a{m,i} rti_qu2{s,m,i}]; 
          Ft_mi{m,i} = [Ft_mi{m,i}; ...
            table(s*ones(size(Ft_smi{s,m,i},1),1), ...
            rti_qu1{s,m,i}', rti_qu2{s,m,i}', ...
            'VariableNames', {'sid', 'rti1', 'rti2'}), ...
            Ft_smi{s,m,i}];
          Fma{m,i} = [Fma{m,i}; Fm{s,m,i}]; 
        elseif s > 6  % send E_37 to J Vision ?                                       
          for nz = 1:2
            rti_qu1a{m,i} = [rti_qu1a{m,i} rti_qu1{s,m,i,nz}]; 
            rti_qu2a{m,i} = [rti_qu2a{m,i} rti_qu2{s,m,i,nz}]; 
            Ft_mi{m,i} = [Ft_mi{m,i}; ...
              table(s*ones(size(Ft_smi{s,m,i,nz},1), 1), ...
              rti_qu1{s,m,i,nz}', rti_qu2{s,m,i,nz}', ...
              'VariableNames', {'sid', 'rti1', 'rti2'}), ...
              Ft_smi{s,m,i,nz}];
            Fma{m,i} = [Fma{m,i}; Fm{s,m,i,nz}]; 
          end
        end                  
      end % s           
    end % m
  end % i 
end 