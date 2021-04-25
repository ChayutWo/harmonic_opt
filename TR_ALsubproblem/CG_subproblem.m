function [output] = CG_subproblem(x,l,u, G, c)
%Use CG to solve subproblem
%min q(x) = 1/2xTGx+xTc st Ax = b
%get current active set
Ax =  getActiveSet(x,l,u);
if length(Ax) == length(x)
    %all constraints are active. Don't have to solve this step
    output = x;
    return;
end   
%Construct matrix A for equality constraint
if isempty(Ax)
    %unconstrain optimization. Solve with normal CG
    P = eye(length(x), length(x));
else
    %equality constrained optimization
    A = zeros(length(Ax),length(x));
    for row = 1:length(Ax)
        A(row, Ax(row)) = 1;
    end
    %construction projection matrix P
    P = eye(length(x), length(x)) - A'*((A*A')\A);
end

% x is already feasible
r = G*x+c;
g = P*r;
d = -g;
tol = 1e-6;
max_iter = 1000;
count = 1;
output = x;
while abs(r'*g) > tol && count < max_iter
    if d'*G*d <=0
        %negative curvature is found
        %move along d until hit some boundary
        [t_bound,~] = calculate_t_bound(x,l,u,-d);
        t_bound = sort(t_bound(t_bound~=0));
        t_min = t_bound(1);
        output = project(x+t_min*d,l,u);
        return
    end
    alpha = r'*g/(d'*G*d);
    [t_bound,~] = calculate_t_bound(x,l,u,-d);
    t_bound = sort(t_bound(t_bound~=0));
    t_min = t_bound(1);
    if alpha > t_min
        %if the update would hit the boundary, stop at the boundary
        output = project(x+t_min*d,l,u);
        return
    end
    x = x + alpha*d;
    r_new = r+alpha*G*d;
    g_new = P*r_new;
    beta = r_new'*g_new/(r'*g);
    d = -g_new+beta*d;
    g = g_new;
    r = r_new;
    count = count + 1;
end
end

