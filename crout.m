function x = crout(a, c, d, b)
% This function solves a linear system Ax = b using Crout's decomposition 
% method for a tridiagonal matrix A.
% The input vectors a, c, and d represent the diagonal, subdiagonal, and 
% superdiagonal elements of the tridiagonal matrix A, respectively.
% Vector b represents the right-hand side of the linear system.

    n = length(a);  % Length of the main diagonal vector a
    n1 = length(c); % Length of the subdiagonal vector c
    n2 = length(d); % Length of the superdiagonal vector d

    % Error checking
    if n1 ~= n2 % The vectors c and d must have the same length
        error('MATLAB:Crout: Not a tridiagonal matrix, incorrect number of elements in parameter arrays.');
    elseif n ~= n1 + 1
        error('MATLAB:Crout: Not a tridiagonal matrix, incorrect number of elements in parameter arrays.');
    end
   
    % Initialize vectors
    p = 1:n;   % Vector to store modified diagonal elements
    q = 1:n-1; % Vector to store modified superdiagonal elements
    x = 1:n;   % Solution vector x
    y = 1:n;   % Intermediate solution vector y
   
    % Main algorithm for Crout's decomposition (forward elimination)
    p(1) = a(1);
    for i = 1:n-1
        q(i) = c(i) / p(i);
        p(i+1) = a(i+1) - d(i) * q(i); % Update diagonal elements, d index adjusted from 1 to n-1
    end
    
    % Forward substitution to solve for y
    y(1) = b(1) / p(1); % y is stored in x for the solution process
    for i = 2:n
        y(i) = (b(i) - d(i-1) * y(i-1)) / p(i);
    end
    
    % Backward substitution to solve for x
    x(n) = y(n);
    for i = (n-1):-1:1
        x(i) = y(i) - q(i) * x(i+1);
    end
end
