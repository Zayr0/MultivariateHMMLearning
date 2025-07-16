% function Kr = krushkal_rank(M)
%     R = rank(M);
%     Kr = R;
% 
%     for i = 1:size(M, 2)
%         v = M(:, i); % Column vector to check
%         S = M(:, 1:end ~= i); % Matrix M with column vector excluded
% 
%         A = sym("a",[size(S, 2) 1]); % Symbolic variables
%         eqns = (sum(S*A)-v == 0); % Create symbolic equation to solving
%         b = struct2cell(solve(eqns,A, Real=true)); % Solve equation
%         n = length(find([b{:}]~=0)); % Find the number of vectors that v is linearly dependent on
% 
%         if n < Kr && n ~= 0
%             Kr = n; % Kruskal rank is equal to the minimum n, for all linearly dependent vectors v_i
%         end
%     end
% end

function k = krushkal_rank(A)
    % Computes the Kruskal rank of a given matrix A
    % The Kruskal rank is the maximum number r such that every subset of r columns is linearly independent.
    
    [~, n] = size(A);
    k = 0;
    
    for r = 1:n
        combs = nchoosek(1:n, r); % All column combinations of size r
        all_independent = true;
        
        for i = 1:size(combs, 1)
            submatrix = A(:, combs(i, :));
            if rank(submatrix) < r
                all_independent = false;
                break;
            end
        end
        
        if all_independent
            k = r;
        else
            break;
        end
    end
end
