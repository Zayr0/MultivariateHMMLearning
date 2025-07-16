function [C,U1,U2,U3]=truncated_mlsvd_3d(T, k)
T1 = mode_n_matricization(T, 1);
T2 = mode_n_matricization(T, 2);
T3 = mode_n_matricization(T, 3);


[U1,~,~] = svd(T1,"econ");
[U2,~,~] = svd(T2,"econ");
[U3,~,~] = svd(T3,"econ");

U1 = U1(:, 1:k);
U2 = U2(:, 1:k);
U3 = U3(:, 1:k);

%% PART 2: Compute the core tensor

% % YOUR CODE GOES HERE
% display(size(T))
% display(size(U1))
% display(size(U2))
% display(size(U3))
C = mode_n_product(T, U1', 1);
C = mode_n_product(C, U2', 2);
C = mode_n_product(C, U3', 3);

end