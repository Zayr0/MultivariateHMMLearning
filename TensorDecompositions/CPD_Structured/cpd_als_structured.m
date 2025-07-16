function [A, B, C, output] = cpd_als_structured(P, k, T_init, O_init, options)

    [d, ~, ~] = size(P);
    
    pi_init = findStationaryDistribution(T_init')';
    
    A = O_init * diag(pi_init) * T_init' * inv(diag(T_init*pi_init));
    B = O_init;
    C = O_init * T_init;


    P1 = mode_n_matricization(P, 1);
    P2 = mode_n_matricization(P, 2);
    P3 = mode_n_matricization(P, 3);


    % emissions = exp(-0.5 * ((1:d) - mean_val).^2 / std_val^2);

    output.relerr = ones(1, options.maxIter);


    kr = khatri_rao((1:d)', ones(1, k));
    one = ones(d, 1);

    for iter = 1:options.maxIter
        B = P2 * khatri_rao(C, A) * pinv((C'*C).*(A'*A));
        B = Stochasticize(B);
        
        C = P3 * khatri_rao(B, A) * pinv((B'*B).*(A'*A));
        C = Stochasticize(C);

        A = P1 * khatri_rao(C, B) * pinv((C'*C).*(B'*B));
        A = Stochasticize(A);

        if(iter > options.maxIter - options.distIter)
            if(options.distribution == "Normal")
                [muA, sigmaA] = fit_normal_distribution(A);
                [muB, sigmaB] = fit_normal_distribution(B);
                [muC, sigmaC] = fit_normal_distribution(C);

                for i = 1:k
                    A(:, i) = pdf(options.distribution, 1:d, muA(i), sigmaA(i));
                    B(:, i) = pdf(options.distribution, 1:d, muB(i), sigmaB(i));
                    C(:, i) = pdf(options.distribution, 1:d, muC(i), sigmaC(i));
                end
        
            elseif(options.distribution == "Beta")
                [a1, a2] = fit_beta_distribution(A);
                [b1, b2] = fit_beta_distribution(B);
                [c1, c2] = fit_beta_distribution(C);

                d_scaled = (1:d+1) / (d+1);

                for i = 1:k
                    A(:, i) = pdf(options.distribution, d_scaled(1:end-1), a1(i), a2(i));
                    B(:, i) = pdf(options.distribution, d_scaled(1:end-1), b1(i), b2(i));
                    C(:, i) = pdf(options.distribution, d_scaled(1:end-1), c1(i), c2(i));
                end

            elseif(options.distribution == "Gaussian")
                sigma = 1;

                A = smoothColumns(A, sigma);
                B = smoothColumns(B, sigma);
                C = smoothColumns(C, sigma);
            end
            A = Stochasticize(A);
            B = Stochasticize(B);
            C = Stochasticize(C);
        end

        output.relerr(iter) = norm(P - cpdgen({A, B, C}), "fro");
    end

    output.relerr = output.relerr(1:iter);
end

function SmoothedMatrix = smoothColumns(Matrix, sigma)
    % smoothColumns applies Gaussian smoothing to each column of the input matrix
    % independently, ensuring that rows do not interfere with each other.
    %
    % Inputs:
    %   - Matrix: A d x k matrix with probability density functions in columns.
    %   - sigma: Standard deviation for Gaussian smoothing.
    %
    % Output:
    %   - SmoothedMatrix: The smoothed d x k matrix.
    %
    
    [d, k] = size(Matrix);
    SmoothedMatrix = zeros(d, k);
    
    % Define a 1D Gaussian filter
    windowSize = ceil(3 * sigma); % Use 3 standard deviations for window size
    x = -windowSize:windowSize;
    gaussFilter = exp(-x.^2 / (2 * sigma^2));
    gaussFilter = gaussFilter / sum(gaussFilter); % Normalize filter
    
    % Apply convolution to each column separately
    for col = 1:k
        SmoothedMatrix(:, col) = conv(Matrix(:, col), gaussFilter, 'same');
    end
end
