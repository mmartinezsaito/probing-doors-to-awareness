 function f = getfitfields(s, fn, lvs)
 f = [];
 for l = lvs  
   f = [f getfield(s, {l}, fn)];
 end
 