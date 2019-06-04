%% To run this script, you first need to install CVX    http://cvxr.com/cvx/download/

%% QUESTION: how should we compare (multiple) graphs?

%% First let us consider two graphs. These could be directed, undirect, weighted, non-weighted, etc.
% in this case we generate one bernoulli graph, and another graph which is
% just like the first one with the nodes permuted.
% it is easy to see, when n is small, and with our eyes, that, although the drawings of the graphs might look different,
% the two graphs are actually the same 

n = 5;
A = (rand(n) > 0.5); A = triu(A,1) + triu(A,1)';
random_permutation = randperm(n);
B = A(random_permutation,random_permutation);

subplot(1,2,1);
plot(graph(A));
subplot(1,2,2);
plot(graph(B));

%% however, when n is large, this is not so easy anymore 

n = 10;
A = (rand(n) > 0.3); A = triu(A,1) + triu(A,1)';
random_permutation = randperm(n);
B = A(random_permutation,random_permutation);

subplot(1,2,1);
plot(graph(A));
subplot(1,2,2);
plot(graph(B));

%% if the graphs are small, it is possible to find if it is the case or not, that the graphs are actuall the same, up to a relabling of the nodes.
% This is called the isomorphism problem. It has been recently claimed by
% László Babai to be solvable in quasi-polynomial.
% note that there might be multiple permutations that generate B from A, some different from one actually choosen
n = 6;
A = (rand(n) > 0.5); A = triu(A,1) + triu(A,1)';
random_permutation = randperm(n);
B = A(random_permutation,random_permutation);

disp(random_permutation);
all_perms = perms(1:n);
for i = 1:size(all_perms, 1)
    if (norm( A(all_perms(i,:),all_perms(i,:)) - B,'fro') == 0)
        disp(all_perms(i,:));
    end
end

%% note that if the graphs are not isomorphic, then we will not be able to find a permutation of the nodes of A such that it equals B
n = 6;
A = (rand(n) > 0.5); A = triu(A,1) + triu(A,1)';
B = (rand(n) > 0.5); B = triu(B,1) + triu(B,1)';

all_perms = perms(1:n);
for i = 1:size(all_perms, 1)
    if (norm( A(all_perms(i,:),all_perms(i,:)) - B,'fro') == 0)
        disp(all_perms(i,:));
    end
end

%% the problem of determining if two graphs are isomorphic is not known to be in P or NP
% a more compicated problem, known to be hard, is that of determining how close two graphs are of being isomorphic
% in other words, the problem of finding some cost function D(graph_1,graph_2) such that if D is small, the graphs are more similar, and if D
% is large, the graphs are more dissimilar. The problem is hard for certain
% natural choices of the function D, not for all choices of D.


% there are two classical approachces for comparingn two graph. Alignment based
% methods, and alignment-free methods. As the name implies, alignment based
% methods first try to align the nodes of the two graphs as well as
% possible, via some permutation of the labels of the nodes. Then, once the nodes are
% aligned, it compares the edges between each pair of nodes in the two
% graphs. These methods try to aling the nodes such that, after the alignemnt, edges in the two graphs overlap as much as possible.
% Alignment free methods, do not do this.

%% canonical example of alignment-free method: spectrum

n = 6;
A = (rand(n) > 0.5); A = triu(A,1) + triu(A,1)';
random_permutation = randperm(n);
B = A(random_permutation,random_permutation);
C = (rand(n) > 0.5); C = triu(C,1) + triu(C,1)';

norm(sort(eigs(A)) - sort(eigs(B))) % spectrum difference in isomorphic graphs
norm(sort(eigs(A)) - sort(eigs(C))) % spectrum difference in non-isomorphic graphs

%% note however that it is possible that two graphs are not-isomorphic, and
% however have the same specturm. So, when comparing spectra, we should keep in
% mind what co-spectrality actually implies.
% We have isomorphism => co-spectrality, but not the other way around

start_list_A = [1,2,3,4,5,6,7,8, 9, 11,12, 13 , 14,15,16,17 ];
end_list_A = [2,3,4,5,6,7,8,9,10, 3, 11, 4,    6,14,15,9 ];

start_list_B = [1,2,3,4,5,6,7,8, 9, 11,12, 13 , 14,15,16,17 ];
end_list_B = [2,3,4,5,6,7,8,9,10, 3,  11,   6,     13,  14,  7,9 ];

A = zeros(17);
B = zeros(17);

for i = 1:length(end_list_A)
    A(start_list_A(i),end_list_A(i)) = 1;
    B(start_list_B(i),end_list_B(i)) = 1;
end

A = 0.0 + (A + A' > 0);
B = 0.0 + (B + B' > 0);
subplot(1,2,1);
plot(graph(A))
subplot(1,2,2);
plot(graph(B))

norm(sort(eigs(A)) - sort(eigs(B)))

%% canonical example of alignemnt-based method: minimum edge difference
% Note that if P is a permutation, for example,  [6     4     5     1     3 2], 
% and A is an adjacancy matrix, then doing A(P,P), is the same as doing
% M A M', where M is a permutation matrix build by putting a 1 on the i-th
% row in the P(i)-th column

n = 6;
A = (rand(n) > 0.5); A = triu(A,1) + triu(A,1)';
random_permutation = randperm(n);

M = zeros(n);
for i = 1:n
    M(i, random_permutation(i)) = 1;
end

A(random_permutation,random_permutation) - M*A*M'

%% once two graphs are aligned, comparing their adjacancy matrices is equal to counting the edge differences between the graphs
% hence, we can estimate how similar two grapphs are by evaluating the
% quantity ||A_aligned - B  || = ||A(P,P) - B || = || M A M' -  B ||, where
% P is the permutation that alignts the two graphs and M is the corresponding permutation
% matrix. This quantity SQUARED, gives us TWICE the nummber of edge
% differences between the two graphs.
% Now notice that || M A M' -  B || = || (  M A  -  B (M')^(-1)   )M'  ||
% and, for many matrix norms, permuting the rows or columns does not change
% its value. Hence, || M A M' -  B || = || (  M A  -  B (M')^(-1)   )M'  ||
% = || (  M A  -  B (M')^(-1)   ) || = || B M -  M A|| , since for a
% permutation (M')^(-1) = M
% Therefore, one canonical alignment-based graph comparison method, is to
% solve the following problem    min_{M in permutation matrices} || B M - M A ||

n = 6;
A = (rand(n) > 0.5); A = triu(A,1) + triu(A,1)';
B = (rand(n) > 0.5); B = triu(B,1) + triu(B,1)';

min_val = inf;
min_perm = nan;
all_perms = perms(1:n);


for i = 1:size(all_perms, 1)
    M = zeros(n);
    for j = 1:n
        M(j, all_perms(i,j)) = 1;
    end

    if (  norm( B*M - M*A,'fro')   <  min_val )
        min_val = norm( B*M - M*A,'fro') ;
        min_perm = M;
    end
        
end

disp(0.5 * min_val^2)

subplot(1,3,1);
plot(graph(A))
subplot(1,3,2);
plot(graph(min_perm*A*min_perm'))
subplot(1,3,3);
plot(graph(B))

%% the problem of the above approach is that it is again hard to compute the optimization
% ONE good possible solution is to relax the constraint in     min_{M in permutation matrices} || B M - M A ||
% such that we work with a more tractable set of matrices.
% one such possibility is replacing M with the set of doubly stochastic
% matrices. These are non-negative matrices whose rows and columns sum to
% 1. In other words, we solve  min_{M in doubly stochastic matrices} || B M - M A ||

n = 6;
A = (rand(n) > 0.5); A = triu(A,1) + triu(A,1)';
B = (rand(n) > 0.5); B = triu(B,1) + triu(B,1)';

cvx_begin quiet

    variable M_relaxed(n,n)
    
    minimize   (     norm(B*M_relaxed - M_relaxed*A,'fro')   )
        subject to
            M_relaxed >= 0
            sum(M_relaxed,1) == 1
            sum(M_relaxed,2) == 1

cvx_end

min_val = inf;
min_perm = nan;
all_perms = perms(1:n);

for i = 1:size(all_perms, 1)
    M = zeros(n);
    for j = 1:n
        M(j, all_perms(i,j)) = 1;
    end

    if (  norm( B*M - M*A,'fro')   <  min_val )
        min_val = norm( B*M - M*A,'fro') ;
        min_perm = M;
    end
        
end

%% notice, however, that we do not get permutation matrices as the output.
% some times, an (almost optimal) permutation can be recovered from the relaxed M. 
% Other times, it is not easy to recover a good permutation for the optimal doubly stochastic matrix
% To recover such permutation, we can solve a simple linear program
% (corresponding to a maximum maching problem on a bipartite graph), which always outputs a
% permutaion matrix


cvx_begin quiet

    variable M_relaxed_bacK_to_perm(n,n)
    
    minimize   (   -trace(M_relaxed_bacK_to_perm*(M_relaxed' + 0.00001*rand(n)))     )
    subject to
        subject to
            M_relaxed_bacK_to_perm >= 0
            sum(M_relaxed_bacK_to_perm) == 1
            sum(M_relaxed_bacK_to_perm') == 1
    
cvx_end

subplot(1,3,1);
imagesc(min_perm);
subplot(1,3,2);
imagesc(M_relaxed);
subplot(1,3,3);
imagesc(M_relaxed_bacK_to_perm);

%% another possible relaxation is using orthogonal matrices for the matrices M
% this leads to   min_{M in orthogonal matrices} || B M - M A ||
% in this case, and when A and B are real and symmetric matrices, it turns
% out that the solution can be computed from a simple Singular Value
% Decomposition

n = 6;
A = (rand(n) > 0.5); A = triu(A,1) + triu(A,1)';
B = (rand(n) > 0.5); B = triu(B,1) + triu(B,1)';

[U_A, D_A, V_A ] = svd(A);  % A = U_A D_A V_A'
[U_B, D_B, V_B ] = svd(B);  % B = U_B D_B V_B'

M_relaxed = U_A*U_B';
norm(D_A - D_B,'fro')

% brute force
min_val = inf;
min_perm = nan;
all_perms = perms(1:n);

for i = 1:size(all_perms, 1)
    M = zeros(n);
    for j = 1:n
        M(j, all_perms(i,j)) = 1;
    end

    if (  norm( B*M - M*A,'fro')   <  min_val )
        min_val = norm( B*M - M*A,'fro') ;
        min_perm = M;
    end
        
end

% project onto permutations
cvx_begin quiet

    variable M_relaxed_bacK_to_perm(n,n)
    
    minimize   (   -trace(M_relaxed_bacK_to_perm*(M_relaxed' + 0.00001*rand(n)))     )
    subject to
        subject to
            M_relaxed_bacK_to_perm >= 0
            sum(M_relaxed_bacK_to_perm) == 1
            sum(M_relaxed_bacK_to_perm') == 1
    
cvx_end



subplot(1,3,1);
imagesc(M_relaxed);
subplot(1,3,2);
imagesc(min_perm);
subplot(1,3,3);
imagesc(M_relaxed_bacK_to_perm);

% in this case it is a bit harder to get a good permutation from the optimal orthogonal
% matrix


%% notice that, a very important property that we want when computing D(graph_1,graph_2)
% is that the distance score is a metric in the mathematical sense. In
% other words, we want it to satisfy D >= 0, D = 0 if and only if graph_1 ~ graph_2,
% D(graph_1,graph_2) = D(graph_2,graph_1), and D(graph_1,graph_2)
% + D(graph_2,graph_3) >= D(graph_1,graph_3)
%
% the reason for wanting this is kind of properties is obvious: if graph 1 and graph 2
% are similar, and graph 2 and graph 3 are also similar, then graph 1 and
% graph 3 should also be similar
%
%
% Futhermore, the metric property allows us to speedup several machine
% learning algorithms, and allows us to improve the performance in other
% learning algorithms. See J. Bento, S. Ioannidis, A Family of tractable graph distances, SDM 2018
% and also see S. Safavi, J. Bento, Tractable n-metrics for multiple graphs, ICML 2019 
%
% Both  min_{M in orthogonal matrices} || B M - M A ||   and     min_{M in doubly stochastic matrices} || B M - M A ||
% result in socres that satisfy the metric properties.
% In fact, for many choices of domain for M, this is the case. See the papers
% above
%

A = (rand(n) > 0.5); A = triu(A,1) + triu(A,1)';
B = (rand(n) > 0.5); B = triu(B,1) + triu(B,1)';
C = (rand(n) > 0.5); C = triu(C,1) + triu(C,1)';

cvx_begin quiet
    variable M_relaxed(n,n)
    minimize   (     norm(B*M_relaxed - M_relaxed*A,'fro')   )
        subject to
            M_relaxed >= 0
            sum(M_relaxed) == 1
            sum(M_relaxed') == 1
cvx_end
d_AB = cvx_optval;
cvx_begin quiet
    variable M_relaxed(n,n)
    minimize   (     norm(A*M_relaxed - M_relaxed*C,'fro')   )
        subject to
            M_relaxed >= 0
            sum(M_relaxed) == 1
            sum(M_relaxed') == 1
cvx_end
d_AC = cvx_optval;
cvx_begin quiet 
    variable M_relaxed(n,n)
    minimize   (     norm(B*M_relaxed - M_relaxed*C,'fro')   )
        subject to
            M_relaxed >= 0
            sum(M_relaxed) == 1
            sum(M_relaxed') == 1
cvx_end
d_BC = cvx_optval;

disp(d_AB + d_BC - d_AC); % we can see that the trianle inequality holds. We could also check that all the other metric properties hold

%% another important question is how to compare, or align, multiple graphs simultaneously
%
%
% One naive way of doing it for N graphs G_1,...,G_N is to define
% D(G_1,..., G_N) = sum_{i,j} D(G_i, G_j),
%
% where each D(G_i, G_j) produces a similarity score, possibily with an
% alignment, between graphs i and j
%
% However, when aligning graphs independently two by two, the resulting
% alignemnts might create problems. For example, D(G_1,G_2) might match
% node x in G_1 to y in G_2. And D(G_2,G_3) might match node y in G_2 to node z
% in G_3. But D(G_1,G_3) might not match x in G_1 to z in G_3.
%
% Therefore, when comparing multiple graphs via alignment-based methods, we
% want to make sure that these alignemnts are consistent.
%
%
% Let M^{ij} be the permutation matrix that aligns G_i to G_j. What we want
% in general is that M^{ij} M^{jk} = M^{ik}, and that M^{ii}  = Identity matrix. This is called ** Alignment Consistency **
%
% This extra constrains makes the problem of computing multigraph similarity scores even harder 
%
%
% One way to solve this problem is directly as 
%
% D(G_1,..., G_N) =  min_{M^{ij} in permutations;} sum_{i,j} ||  A_i M^{ij}   - M^{ij} A_j || 
% subjec to M^{ij} M^{jk} = M^{ik}, and that M^{ii}  = Identity matrix.
%
%
% We could then try to relax the constraint that M^{ij} are permutations. 
% D(G_1,..., G_N) =  min_{M^{ij} in doubly stochastic matrices} sum_{i,j} ||  A_i M^{ij}   - M^{ij} A_j || 
% subjec to M^{ij} M^{jk} = M^{ik}, and that M^{ii}  = Identity matrix.
% but the Alignment Consistency makes the problem non convex, and hence hard to solve
%
% Another way to solve this problem is to do compute a distance as follows
%
%
% D(G_1,..., G_N) =  min_{G_ref} D(G_i,G_ref) 
% and this can be further specified as
% D(G_1,..., G_N) =  min_{G_ref} D(G_i,G_ref) =   min_{Q^{i} in permutations; A_ref in adjacancy matrices} sum_{i} ||  A_i Q^{i}   - Q^{i}A_ref || 
%
% This is called the Fermat distance induced by D(.,.).
% This has the advantage that, if we define M^{ij} to be Q^i * inv(Q^j) = Q^i * Q'^j ),
% the resulting M^{ij} are, by definition, alignment consistent. For
% example, M^{ij} M^{jk} = Q^i inv(Q^j) Q^j inv(Q^k) = Q^i inv(Q^k) = M^{ik} 
%
%
% It also has the advantage of directly leading to continuous optimization
% problems, if the constraint that Q^i are permutations are relaxed.
%
%
% Finally, it has the advantage that D(G_1,..., G_N) satisfies some form of
% generalized metric property, including a generalization of the triangle inequality for multiple elements.
%
%
% To obtain a relazation we could, for example, define
% D(G_1,..., G_N) =   min_{Q^{i} in doubly stochastic matrices; A_ref in real matrices with entires in [0,1] } sum_{i,j} ||  A_i Q^{i}   - Q^{i}A_ref || 
%
%
% The probem is that this optimization problem is non-convex, and hence hard to solve in general

%% A different approach is to relax both the permutations, and the alignment consistency constraint
% It is not hard to prove that the constraint  M^{ij} M^{jk} = M^{ik} and
% M^{ii}  = Identity matrix, imply that M^{ij} can be written as Q^i Q^j'
% for some permutation matrices Q^i. This means that, if we define BigM to
% be a block matrix with (ij)-th block equal to M^{ij}, we will have that
% BigM is low rank and positive semi-definite.
% We use these facts to define two multi-graph scores that can be computed
% via convex optimization
%
%
% Here we show how they look for D(G_1, G_2, G_3), comparing k = 3 graphs

    nuc_vs_psd = 1; % flag is 1 if we are imposing the low rank constraint via nuclear norm. And flag is 0, if we are just imposing semi-definiteness
    
    n = 6;
    k = 3;
    
    A = (rand(n) > 0.5); A = triu(A,1) + triu(A,1)';
    B = (rand(n) > 0.5); B = triu(B,1) + triu(B,1)';
    C = (rand(n) > 0.5); C = triu(C,1) + triu(C,1)';

    all_graphs.adj{1} = A;
    all_graphs.adj{2} = B;
    all_graphs.adj{3} = C;
    

    A = zeros(n,n,k);
    for i = 1:k
        A(:,:,i) = all_graphs.adj{i};
    end

    cvx_begin quiet

    variable P(n*k,n*k)  % this define BigM, that we will optimize. A block matrix. Each block being M^{ij} above

    % this defines the objective   sum_{i,j} ||  A_i M^{ij}   - M^{ij} A_j || 
    s = 0;
    for i = 1:k
        for j = 1:k
            s = s + norm(A(:,:,i)*P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) ) - P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) )*A(:,:,j) , 'fro');
        end
    end

    minimize (0.5*s   )

    subject to  
            if (nuc_vs_psd == 1)
                norm_nuc(P) <= n*k;
            else
            	P == semidefinite(n*k);    
            end
            
            diag(P) == 1;
            for i = 1:k
                for j = 1:k
                    P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) ) >= 0;
                    sum(P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) )) == 1;
                    sum(P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) )') == 1;
                end
            end
    cvx_end

    
    imagesc(P)
    
% we can then extract the different matrices M^ij , and we if want, project
% them onto permutations to get 0-1. Note, however, that this projection
% often distorts the alignment consistency property

