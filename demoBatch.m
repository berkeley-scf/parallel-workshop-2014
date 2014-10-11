feature('numThreads', 1); 
NSLOTS = str2num(getenv('NSLOTS')); # if running on the cluster
% NSLOTS = 8; # otherwise choose how many cores you want to use
pool = parpool(NSLOTS); 
% assume you have test.m with a function, test, taking two inputs (n and seed)
% and returning 1 output
n = 10000000;
job = cell(1,6); 
job{1} = parfeval(pool, @test, 1, n, 1);  
job{2} = parfeval(pool, @test, 1, n, 2);  
job{3} = parfeval(pool, @test, 1, n, 3);  
job{4} = parfeval(pool, @test, 1, n, 4);  
job{5} = parfeval(pool, @test, 1, n, 5);  
job{6} = parfeval(pool, @test, 1, n, 6);  

% wait for outputs, in order
output = cell(1, 6);
for idx = 1:6
  output{idx} = fetchOutputs(job{idx});
end 

% alternative way to loop over jobs:
for idx = 1:6
  jobs(idx) = parfeval(pool, @test, 1, n, idx); 
end 

% wait for outputs as they finish
output = cell(1, 6);
for idx = 1:6
  [completedIdx, value] = fetchNext(jobs);
  output{completedIdx} = value;
end 

delete(pool);


%%% w/ threading

NSLOTS = str2num(getenv('NSLOTS')); # if running on the cluster
% NSLOTS = 8; # otherwise choose how many cores you want to use
pool = parpool(NSLOTS); 
n = 5000;
nJobs = 2;
pool = parpool(nJobs);
% pass number of threads as number of slots divided by number of jobs
% testThread() function should then do: 
% feature('numThreads', nThreads);
% where nThreads is the name of the relevant argument to testThread()
jobt1 = parfeval(pool, @testThread, 1, n, 1, NSLOTS/nJobs);
jobt2 = parfeval(pool, @testThread, 1, n, 1, NSLOTS/nJobs);
jobt3 = parfeval(pool, @testThread, 1, n, 1, NSLOTS/nJobs);
jobt4 = parfeval(pool, @testThread, 1, n, 1, NSLOTS/nJobs);

output1 = fetchOutputs(jobt1);
output2 = fetchOutputs(jobt2);
output3 = fetchOutputs(jobt3);
output4 = fetchOutputs(jobt4);

delete(pool);