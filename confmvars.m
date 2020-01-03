function confmvars(cm)

tp = cm(1,1); % hits 
fp = cm(1,2); % false alarm, type I error 
fn = cm(2,1); % miss, type II error  
tn = cm(2,2); % correct rejection
cp = sum(cm(:,1));
cn = sum(cm(:,2));
dp = sum(cm(1,:));
dn = sum(cm(2,:));

sensitivity = tp / (tp + fn)
specificity = tn / (tn + fp)
precision   = tp / (tp + fp)
negpredval  = tn / (tn + fn) % for type II, this equals SDI
accuracy    = (tp + tn) / sum(cm(:))

