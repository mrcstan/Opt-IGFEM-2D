global G_var
G_var = 0;

nCompleted = 0;
inputs = [pi/2,pi];

p = gcp();

%
for i = 1:p.NumWorkers
    futures(i)= parfeval(p,@test_func,1,inputs(i));
end

%out = fetchOutputs(futures);
cancelFutures = onCleanup(@() cancel(futures));
out = cell(p.NumWorkers,1);
while nCompleted < p.NumWorkers
    [completedIdx, resultThisTime] = fetchNext(futures);
    if (~isempty(completedIdx))
        
        nCompleted = nCompleted + 1;
        out{completedIdx} = resultThisTime;
    end
end
%
%{
futures = parfevalOnAll(p,@test_func,1,inputs(1));

out = fetchOutputs(futures);
%}