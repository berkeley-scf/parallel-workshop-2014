nslots = str2num(getenv('NSLOTS')); % if running on the cluster
% nslots = 8; # otherwise choose how many cores you want to use
mypool = parpool(nslots); 
% parpool open local nslots # alternative

n = 3000
nIts = 500
c = zeros(n, nIts);
parfor i = 1:nIts
     c(:,i) = eig(rand(n)); 
end

delete(mypool) 
% delete(gcp) works if you forget to name your pool by assigning the 
% output of parpool to a variable

