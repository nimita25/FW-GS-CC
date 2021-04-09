function [fn_val, sol_SDP] = CC_SDPT(filename, e1, MemLog)
% CC_SDPT(filename, epsilon, MemLog)
% where filename is the name of the '.mat' file containing information of
% input graph,
% MemLog = 1 if you want to track the memory usage
% CC_SDPT solves the problem to epsilon-optimality

tic;
%rng(10);

%% Check input data
narginchk(2,3);
if isempty(MemLog), MemLog = 0; end

%% Track memory usage
warning off all;
profile clear
if MemLog == 1, profile -memory on; end

%% Generate input graph based on input type
if isfile(['CCGraphs/',filename])
    [V,L,W] = CreateInputGraph(filename);
else
    disp('Enter valid filename');
    return
end
E = CreateEdgeSet(V,L,W);
%Cost matrix of SDP relaxation
C = L+W;

%% Set input parameters
n = V;
n2 = length(E);
e = ones(n,1); 
C1{1} = -C; 
b = e;
blk{1,1} = 's';  blk{1,2} = n;
blk{2,1} = 'l';  blk{2,2} = n2; 

A = cell(1,n+n2);
for i = 1:n; A{i} = sparse(i,i,1,n,n); end
cnt = n+1;
for i = 1:n2
    A{1,cnt} = sparse(E(i,1), E(i,2),1,n,n);
    cnt = cnt+1;
end

%Parameters for inequality constraints
Avec = svec(blk(1,:),A,1);
Avec{2,1} = [sparse(n2,n), -speye(n2,n2)];
b = [b; zeros(n2,1)]; 
C1{2,1} = zeros(n2,1);

%% Solve the problem
OPTIONS.gaptol = e1;
%[obj,sol_SDP,y,Z] = sqlp(blk,Avec,C1,b,OPTIONS,X0,y0,Z0);
[obj,sol_SDP,y,Z] = sqlp(blk,Avec,C1,b,OPTIONS);
fn_val = -obj(1);
toc
time = toc;

%% Compute feasible solution with best value
nsamples = 2;
Z = mvnrnd(zeros(n,1),sol_SDP{1},nsamples);
Z = Z';
disp(size(Z));
%Assign nodes to partition
assignment = zeros(V,1);
for i = 1:V
    if sign(Z(i,1)) == sign(Z(i,2))
        if sign(Z(i,1)) < 0
            assignment(i) = 1;
        else
            assignment(i) = 4;
        end
    else
        if sign(Z(i,1)) < 0
            assignment(i) = 2;
        else
            assignment(i) = 3;
        end
    end
end
%Generate the max-agree objective value for the generated cluster
e = @(k,n)[zeros(k-1,1);1;zeros(n-k,1)];
Cut = 0;
for i = 1:n2
    if assignment(E(i,1)) ~= assignment(E(i,2))
        Cut = Cut - e(E(i,1),V)'*L*e(E(i,2),V);
    else
        Cut = Cut + e(E(i,1),V)'*W*e(E(i,2),V);
    end
end
Cut = Cut*2;


%% Display output
disp(['SDP objective function value:',num2str(fn_val)]);
disp(['Value of the generated cluster:',num2str(Cut)]);
if MemLog == 1
    p = profile('info');
    memoryUsed = max([p.FunctionTable.PeakMem]);
    disp(['Memory used (in kB):',num2str(memoryUsed/1024)]);
end
disp(['Time (in secs):',num2str(time)]);

%% Save output
if ~exist('output-SDPT','dir'), mkdir('output-SDPT'); end
CC.SDPT.SDPObj = fn_val;
CC.SDPT.ClusterValue = Cut;
if MemLog == 1
    p = profile('info');
    memoryUsed = max([p.FunctionTable.PeakMem]);
    memoryUsed = [num2str(memoryUsed/1024),' kB'];
else
    memoryUsed = 'Memory not logged';
end
CC.SDPT.MemoryUsed = memoryUsed;
CC.SDPT.Time = time;

save(['output-SDPT/',filename],'CC','-v7.3');

%% Functions used in the file
    function [V,L,W] = CreateInputGraph(InputData)
        p = load(['CCGraphs/',InputData]);
        V = size(p.Problem.W1,1);
        L = spdiags(p.Problem.W1*ones(V,1),0,V,V) - p.Problem.W1;
        W = p.Problem.W2;
            
    end


    function [E] = CreateEdgeSet(V,L,W)
        %%%%Generate the list of edges E^+ = (i,j), i<j
        E1 = zeros(nnz(W),2);
        cnt = 1;
        for i = 1:V
            idx = find(W(i,i+1:V));
            idx = idx+i;      %% adjust index
            for j = 1:length(idx)
                E1(cnt,:) = [i,idx(j)];
                cnt = cnt+1;
            end
        end
        E1 = E1(any(E1,2),:);

        %%%%Generate the list of edges E^- = (i,j), i<j
        E2 = zeros(nnz(L),2);
        cnt = 1;
        for i = 1:V
            idx = find(L(i,i+1:V));
            idx = idx+i;      %% adjust index
            for j = 1:length(idx)
                E2(cnt,:) = [i,idx(j)];
                cnt = cnt+1;
            end
        end
        E2 = E2(any(E2,2),:);
        E = [E1;E2];
        E = unique(E,'rows');
    end


end
