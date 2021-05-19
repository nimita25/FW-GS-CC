function [] = CorClus(InputType, InputData, e1, max_time, MemLog)
% CorClus(InputType, InputData, e1, max_time, MemLog)
% (a) If you want to solve CC on a random graph
% CorClus('R', [V,degree], e1, max_time, MemLog)
% where V number of nodes and average degree of each node equal to 'degree'
% 
% (b) If you want to solve CC on the input graph
% CorClus('S',filename, e1, max_time, MemLog)
% where filename is the '.mat' that consists of graph information in the format
% given in Step 1.
% 
% Set e1 = epsilon, the relative error to generate the solution to SDP.
% Optionally, set
% (i) MemLog = 1 if you want to track memory usage
% (ii) max_time = max time to run the algorithm (in seconds)
% The output is stored in the 'output_CC' directory


%% Set seed value
rng(10);
tic;

%% Check input data
narginchk(3,5);
if nargin == 3
    max_time = 36000;
    MemLog = 0;
end
if nargin == 4
    MemLog = 0;
end

%% Track memory usage
warning off all;
profile clear
if MemLog == 1
    profile -memory on;
end

%% Set parameters for computing eigenvector
maxit = 300;
tolerance = 0.1;
opts.isreal = 1;
opts.issym = 1;
opts.maxit = maxit;
opts.tol = tolerance;

%% Generate input graph based on input type
if isequal(InputType,'R')
    if ~isnumeric(InputData(1)) && ~isnumeric(InputData(2))
        disp('Enter valid number of nodes (V) and degree');
        return    
    end
    [V,L,W,W1] = CreateRandomGraph(InputData);
elseif isequal(InputType,'S')
    if isfile(['CCGraphs/',InputData])
        disp('Reading input data');
        [V,L,W] = CreateInputGraph(InputData);
    else
        disp('Enter valid filename');
        return
    end
else
    disp('Enter valid Input type: R or S');
    return
end
[L,W] = ShuffleNodes(V,L,W);
E = CreateEdgeSet(V,L,W);
%Cost matrix of SDP relaxation
C = L+W;

%% Set input parameters for FW-GS
disp('Setting parameters');
alpha = V;
nedges = length(E);
nsamples = 10;
c = 4;
M = (c*log(2*V+nedges))/e1;
beta = c*(trace(L)+sum(sum(W)));
Cf = beta*M*alpha^2;
eta = 10^-4;
lambdamaxC = eigs(C,1,'LR');

%% Initialize the algorithm with identity matrix
%Create matrix of random samples
disp('Creating initial random samples');
x = 1:V;
for ii = 1:2*nsamples
    z = randn(V,1);
    if ii == 1
        Z = z;
    else
        Z = [Z z];
    end
end
%Create linear mapping of objective function [<C,X>, A(X)]
disp('Creating linear mapping');
t = 1;
gamma = 2/(t+2);
v = [ones(V+1,1);zeros(nedges,1)];
v(1) = trace(L);

%% First iteration of the algorithm
disp('Starting FW-GS');
%Compute the gradient of obj function
max_viol = max([abs(M.*(v(2:V+1)-1));M.*(-v(V+2:length(v)))]);
penalty = cat(1,exp(M.*(v(2:V+1)-1)-max_viol),exp(M.*(1-v(2:V+1))-max_viol),exp(M.*(-v(V+2:length(v)))-max_viol));
matrix_entry = cat(1,(beta/sum(penalty))*(penalty(1:V)-penalty(V+1:2*V)), -(beta/sum(penalty))*penalty(2*V+1:length(penalty)));
%Generate the update direction using 'eigs'
[u,l] = eigs(-sparse([x,E(:,1)',E(:,2)'],[x,E(:,2)',E(:,1)'],[matrix_entry(1:V);matrix_entry(V+1:V+nedges)/2;matrix_entry(V+1:V+nedges)/2])+C,1,'LR',opts);
%Compute update direction
h = zeros(V+1+nedges,1);
if l >= 0
    h(1) = alpha*u'*C*u;
    h(2:V+1) = alpha*u.^2;
    for ii = 1:nedges
        h(V+1+ii) = alpha*u(E(ii,1))*u(E(ii,2));
    end
end
%Initial constraint violation
%vecMaxPenalty = [];
sumW = sum(sum(W));

%% Run the loop until gap is less than epsilon*(Tr(C)+sum(W))
disp('#itr|Duality Gap|Max Constr Viol');
while ( ((h-v)'*[1;-matrix_entry] > e1*(trace(C)+sumW) ) && toc <= max_time)
    %Display intermediate status
    if mod(t,1000) == 0
        disp([int2str(t),'|',num2str(round((h-v)'*[1;-matrix_entry],2)),'|',num2str(round(max_viol/M,2))]);
    end
    
    %% Generate samples from gradient and update the samples
    for ii = 1:2*nsamples
        x1 = normrnd(0,1);
        w = u*x1;
        if l >= 0
            Z(:,ii) = sqrt(1-gamma)*Z(:,ii) + sqrt(gamma)*w;
        else
            Z(:,ii) = sqrt(1-gamma)*Z(:,ii);
        end
    end
    
    %% Update the linear map 'v'
    v = (1-gamma)*v + gamma*h;
    t = t+1;
    gamma = 2/(t+2);
    
    
    %% Compute the update direction using 'eigs'
    %Compute gradient at the new point
    max_viol = max([abs(M.*(v(2:V+1)-1));M.*(-v(V+2:length(v)))]);
    penalty = cat(1,exp(M.*(v(2:V+1)-1)-max_viol),exp(M.*(1-v(2:V+1))-max_viol),exp(M.*(-v(V+2:length(v)))-max_viol));
    matrix_entry = cat(1,(beta/sum(penalty))*(penalty(1:V)-penalty(V+1:2*V)), -(beta/sum(penalty))*penalty(2*V+1:length(penalty)));
    %Compute eigenvector using eigs
    %vecMaxPenalty = [vecMaxPenalty; max_viol];
    [u,l] = eigs(-sparse([x,E(:,1)',E(:,2)'],[x,E(:,2)',E(:,1)'],[matrix_entry(1:V);matrix_entry(V+1:V+nedges)/2;matrix_entry(V+1:V+nedges)/2])+C,1,'LR',opts);
    %Compute update direction
    h = zeros(V+1+nedges,1);
    if l >= 0
        h(1) = alpha*u'*C*u;
        h(2:V+1) = alpha*u.^2;
        for ii = 1:nedges
            h(V+1+ii) = alpha*u(E(ii,1))*u(E(ii,2));
        end
    end
    
end

%% Compute feasible solution with best value
%Generate feasible samples
ineq = max(max(-v(V+2:length(v))),0);
maxeq = max(v(2:V+1))+ineq;
for ii = 1:2*nsamples
    y = randn;
    z = Z(:,ii)+sqrt(ineq)*y*ones(V,1);
    rndm = randn(V,1);
    w = sqrt(ones(V,1)-((v(2:V+1)+ineq)/maxeq)).*rndm+z/sqrt(maxeq);
    Z(:,ii) = w;
end
%Assign nodes to partition
assignment = zeros(V,nsamples);
for jj = 1:nsamples
    for ii = 1:V
        if sign(Z(ii,2*(jj-1)+1)) == sign(Z(ii,2*jj))
            if sign(Z(ii,2*(jj-1)+1)) < 0
                assignment(ii,jj) = 1;
            else
                assignment(ii,jj) = 4;
            end
        else
            if sign(Z(ii,2*(jj-1)+1)) < 0
                assignment(ii,jj) = 2;
            else
                assignment(ii,jj) = 3;
            end
        end
    end
end
%Generate the max-agree objective value for the generated cluster
e = @(k,n)[zeros(k-1,1);1;zeros(n-k,1)];
BestCut = 0;
AvgCut = 0;
for jj = 1:nsamples
    Cut = 0;
    for ii = 1:nedges
        if assignment(E(ii,1),jj) ~= assignment(E(ii,2),jj)
            Cut = Cut - e(E(ii,1),V)'*L*e(E(ii,2),V);
        else
            Cut = Cut + e(E(ii,1),V)'*W*e(E(ii,2),V);
        end
    end
    Cut = Cut*2;
    AvgCut = AvgCut+Cut;
    if Cut >= BestCut
        BestCut = Cut;
    end
end
AvgCut = AvgCut/nsamples;
disp(BestCut);

toc
time = toc;

%% Display output
disp('Bound on duality gap (stopping criteria):');
disp((h-v)'*[1;-matrix_entry]);
sp = v(1);
constr_viol_eq = norm(v(2:V+1) -ones(V,1),inf);
constr_viol_ineq = max(-v(V+2:length(v)));

disp('Number of iterations:');
disp(t);

disp('Value of SDP relaxation');
disp(sp);

disp('Maximum violation of constraints (equality,inequality):');
disp(constr_viol_eq);
disp(constr_viol_ineq);

disp('Solution of max-agreement:');
disp(BestCut);
disp(nedges);

%% Write output

CC.InputParams.nNodes = V;
CC.InputParams.nEdges = nedges;
CC.InputParams.epsilon = e1;
CC.InputParams.StopCrit = e1*(trace(C)+sumW);
CC.InputParams.MaxTime = max_time;
CC.Output.SDPObjVal = sp;
CC.Output.MaxInfeasIneq = max(0,constr_viol_ineq);
CC.Output.MaxInfeasEq = constr_viol_eq;
CC.Output.BestClusterValue = BestCut;
CC.Output.AvgClusterValue = AvgCut;
CC.NIterations = t;
CC.Time = time;
if MemLog == 1
    p = profile('info');
    memoryUsed = max([p.FunctionTable.PeakMem]);
    memoryUsed = [num2str(memoryUsed/1024),' kB'];
else
    memoryUsed = 'Memory not logged';
end
CC.MemoryUsed = memoryUsed;
if (h-v)'*[1;-matrix_entry] > e1*(trace(C)+sumW)
    CC.Status = 'Maximum time reached';
else
    CC.Status = 'Approx solution found';
end 

if ~exist('output_CC','dir'), mkdir('output_CC'); end
if InputType == 'R'
    filename = ['R-',datestr(now,'dd-mm-yy-HH:MM-'),int2str(V),'-',int2str(nedges)];
    CC.InputParams.Data = ['Random graph with ',int2str(V),' nodes and ',int2str(nedges),' edges'];
    Problem.W1 = W1;
    Problem.W2 = W;
    save(['CCGraphs/R-',datestr(now,'dd-mm-yy-HH:MM-'),int2str(V),'-',int2str(length(E)),'.mat'],'Problem');
else
    CC.InputParams.Data = InputData;
    filename = [datestr(now,'dd-mm-yy-HH:MM-'),InputData];
end
save(['output_CC/',filename],'CC','-v7.3');

%% Functions used in the file
    function [V,L,W,W1] = CreateRandomGraph(InputData)
        disp('Creating two random graphs');
        
        V = InputData(1);
        %Create random graph of similarity 'W'
        Edges = floor(V*InputData(2)); %Density of similarity graph
        idx = randi(V,Edges,2);
        %Create edges to make the graph connected
        s = floor(rand(1,V-1).*(1:V-1))+1;
        t = 2:V;
        idx1 = [s',t'];
        %idx2 - matrix of all the edges of random connected graph
        idx2 = [idx;idx1];
        clear s
        clear t
        clear idx1
        clear idx
        %Create a symmetric matrix
        idx2 = [idx2; [idx2(:,2),idx2(:,1)]];
        %Delete repeated edges
        idx2 = unique(idx2,'rows');
        %Delete self-loops
        idx2(idx2(:,1)==idx2(:,2),:) = [];
        val = ones(length(idx2),1); %nonnegative weight of similar edges
        %Create a sparse matrx
        W = sparse(idx2(:,1),idx2(:,2), val/2,V,V);

        %Create random graph of dissimilarity 'L'
        Edges = floor(V*InputData(2)); %Density of dissimilarity graph
        idx = randi(V,Edges,2);
        %Create edges to make the graph connected
        s = floor(rand(1,V-1).*(1:V-1))+1;
        t = 2:V;
        idx1 = [s',t'];
        %idx2 - matrix of all the edges of random connected graph
        idx2 = [idx;idx1];
        clear s
        clear t
        clear idx1
        clear idx
        %Create a symmetric matrix
        idx2 = [idx2; [idx2(:,2),idx2(:,1)]];
        %Delete repeated edges
        idx2 = unique(idx2,'rows');
        %Delete self-loops
        idx2(idx2(:,1)==idx2(:,2),:) = [];
        val = ones(length(idx2),1); %nonnegative weight of similar edges
        %Create a sparse matrx
        W1 = sparse(idx2(:,1),idx2(:,2), val/2,V,V);
        clear a
        clear b
        %Create a sparse matrix
        L = spdiags(W1*ones(V,1),0,V,V) - W1;
    end

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

    function [L1,W1] = ShuffleNodes(V,L,W)
        nodes = randperm(V);
        ll = zeros(nnz(L),3);
        cnt = 1;
        for i = 1:V
            idx = find(L(i,i+1:V));
            idx = idx+i;
            for j = 1:length(idx)
                ll(cnt, 1) = nodes(i);
                ll(cnt, 2) = nodes(idx(j));
                ll(cnt, 3) = -L(i,idx(j));
                cnt = cnt+1;
                ll(cnt, 2) = nodes(i);
                ll(cnt, 1) = nodes(idx(j));
                ll(cnt, 3) = -L(i,idx(j));
                cnt = cnt+1;
            end
        end
        
        ll(ll(:,1)==0,:) = [];
        L2 = sparse(ll(:,1),ll(:,2), ll(:,3),V,V);
        L1 = spdiags(L2*ones(V,1),0,V,V) - L2;
        
        
        ll = zeros(nnz(W),3);
        cnt = 1;
        for i = 1:V
            idx = find(W(i,i+1:V));
            idx = idx+i;
            for j = 1:length(idx)
                ll(cnt, 1) = nodes(i);
                ll(cnt, 2) = nodes(idx(j));
                ll(cnt, 3) = W(i,idx(j));
                cnt = cnt+1;
                ll(cnt, 2) = nodes(i);
                ll(cnt, 1) = nodes(idx(j));
                ll(cnt, 3) = W(i,idx(j));
                cnt = cnt+1;
            end
        end
        W1 = sparse(ll(:,1),ll(:,2),ll(:,3),V,V);
        
        
                
    end



end
