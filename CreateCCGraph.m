function [] = CreateCCGraph(filename)

disp('Loading input file');
delta = 0.05;
p = load(['Datasets/',filename]);
V = size(p.Problem.A,1);
C = p.Problem.A;
C = 0.5*(C+C'); %Create symmetric matrix

clear LEdges
for i = 1:V
    idx = find(C(i,1:V));
    LEdges{i} = idx;
end

Lidx1 = zeros(V*(V-1),0);
Lidx2 = zeros(V*(V-1),0);
Lval = zeros(V*(V-1),0);
Widx1 = zeros(V*(V-1),0);
Widx2 = zeros(V*(V-1),0);
Wval = zeros(V*(V-1),0);
Lcnt = 0;
Wcnt = 0;
disp('Creating the two graphs');
for i = 1:V
    for j = i+1:V
        J = numel(intersect(LEdges{i},LEdges{j}))/numel(union(LEdges{i},LEdges{j}));
        S = log((1+J-delta)/(1-J+delta));
        if S > -0.04 && S < 0
            Lcnt = Lcnt + 1;
            Lidx1(Lcnt) = i;
            Lidx2(Lcnt) = j;
            Lval(Lcnt) = -S;
        elseif S > 0.2
            Wcnt = Wcnt + 1;
            Widx1(Wcnt) = i;
            Widx2(Wcnt) = j;
            Wval(Wcnt) = S;
        end
    end
end
L = sparse(Lidx1,Lidx2, Lval/2,V,V);
W = sparse(Widx1,Widx2, Wval/2,V,V);
L = 0.5*(L+L');
W = 0.5*(W+W');

Problem.W1 = L;
Problem.W2 = W;
if ~exist('CCGraphs','dir'), mkdir('CCGraphs'); end
save(['CCGraphs/CC-',filename],'Problem');


end
