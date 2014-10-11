function out = testThread(n, seed, nThreads)
feature('numThreads', nThreads);
% when on the cluster, note you need to ensure that when considering tasks and threads you don't have more than NSLOTS cores occupied at once
rng(seed);

mat = normrnd(0, 1, n, n);

for i = 1:100
  sigma = mat'*mat;
  L = chol(sigma);
end

out = L(1,1);

end