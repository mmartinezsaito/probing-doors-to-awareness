function [facvec, levs] = make_facvec(Fm_smif, s, varargin)
facvec = Fm_smif;
if s < 7
    fvnans = isnan(facvec);
    levs = unique(facvec(~fvnans));
    nl = length(levs);
    facvec(fvnans) = levs(randi(nl));
    for li = 1:nl, ni(li) = sum(facvec == levs(li)); end
else
    facvec1 = varargin{1};
    fi = varargin{2};
    facvecc = {facvec facvec1};
    levs = {}; nl = {};
    for j = 1:2
        fvnans = isnan(facvecc{j});
        levs{j} = unique(facvecc{j}(~fvnans));
        nl{j} = length(levs{j});
        facvecc{j}(fvnans) = levs{j}(randi(nl{j}));
        for li = 1:nl{j}
            ni(li,j) = sum(facvecc{j} == levs{j}(li));
        end
    end
    switch fi
        case 5,    levs = [0 1];
        case 4,    levs = [-4:-1 1:4];
        otherwise, levs = union(levs{1}, levs{2});
    end
    ni = sum(ni, 2);
    nl = length(levs);
    facvec = facvecc;
end