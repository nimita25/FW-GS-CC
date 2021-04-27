function [] = CreateCCGraph_Labels(filename)

disp('Loading input file');
p = load(['Datasets/',filename]);
V = size(p.Problem.A,1);
C = p.Problem.A;
C = 0.5*(C+C'); %Create symmetric matrix

%% Check whether the edges of the input graph are signed
% If negative-weight edges exist, they are labelled are 'dissimilar' and
% the rest as 'similar
% If the edges of the input graph have nonnegative weights, generate labels
% and weights for edge based on Jaccard coefficient
if nnz(C<0) ~= 0
    mask = C<0;
    L = sparse(V,V);
    L(mask) = -C(mask);
    mask  = C>0;
    W = sparse(V,V);
    W(mask) = C(mask);
else
    %% Generate the labels based on Jaccard coefficient
    disp('Creating the two graphs');
    clear LEdges
    LEdges = cell(V,1);
    for i = 1:V
        idx = find(C(i,1:V));
        LEdges{i} = idx;
    end

    delta = 0.05;
    Lidx1 = zeros(nnz(C),0);
    Lidx2 = zeros(nnz(C),0);
    Lval = zeros(nnz(C),0);
    Widx1 = zeros(nnz(C),0);
    Widx2 = zeros(nnz(C),0);
    Wval = zeros(nnz(C),0);
    Lcnt = 0;
    Wcnt = 0;
    for i = 1:V
        for j = LEdges{i}
            J = numel(intersect(LEdges{i},LEdges{j}))/numel(union(LEdges{i},LEdges{j}));
            S = log((1+J-delta)/(1-J+delta));
            if S < 0
                Lcnt = Lcnt + 1;
                Lidx1(Lcnt) = i;
                Lidx2(Lcnt) = j;
                Lval(Lcnt) = -S;
            else
                Wcnt = Wcnt + 1;
                Widx1(Wcnt) = i;
                Widx2(Wcnt) = j;
                Wval(Wcnt) = S;
            end
        end
    end
    L = sparse(Lidx1,Lidx2, Lval/2,V,V);
    W = sparse(Widx1,Widx2, Wval/2,V,V);
end
            

%% Save the file
L = 0.5*(L+L'); %Create symmetric matrix
W = 0.5*(W+W');
disp('Saving file');
Problem.W1 = L;
Problem.W2 = W;
if ~exist('CCGraphs','dir'), mkdir('CCGraphs'); end
save(['CCGraphs/L-',filename],'Problem');


end
