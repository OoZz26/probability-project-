% % Set the parameters
% num_realizations = 30;
% t = 0:0.01:2; % Time vector from 0 to 2 with step 0.1
% 
% % Initialize the matrix to store realizations
% matrix = zeros(num_realizations, length(t));
% 
% % Generate realizations
% for i = 1:num_realizations
%     theta = unifrnd(-pi, pi); % Generate random theta from a uniform distribution
%     Z_t = cos(4 * pi * t + theta);
%     matrix(i, :) = Z_t;
% end
% 
% % Print or use the matrix as needed
% disp(matrix);
% plot(t,matrix(8, :))
% 
% 
% 
% Generate the matrix of realizations
num_realizations = 30;
time_vector = 0:0.01:2;

matrix_Z = zeros(num_realizations, length(time_vector));

for i = 1:num_realizations
    theta = unifrnd(-pi, pi);
    Z_t = cos(4 * pi * time_vector + theta);
    matrix_Z(i, :) = Z_t;
end

% Define autocorrelation function p(xi, xi)
p_autocorrelation = @(xi) exp(-xi^2);

% Evaluate the cross-correlation function
R_x = zeros(length(time_vector));

for ti_idx = 1:length(time_vector)
    for tj_idx = 1:length(time_vector)
        ti = time_vector(ti_idx);
        tj = time_vector(tj_idx);
        
        integrand = zeros(1, num_realizations);
        for k = 1:num_realizations
            integrand(k) = matrix_Z(k, ti_idx) * matrix_Z(k, tj_idx) * p_autocorrelation(matrix_Z(k, ti_idx));
        end
        
        R_x(ti_idx, tj_idx) = trapz(integrand);
    end
end
hold on



% 
% disp('Cross-Correlation Matrix:');
% disp(R_x);
% 
% % plot(time_vector,R_x(1,:))
% % Plot 3D surface
% figure;
% surf(time_vector, time_vector, R_x);
% xlabel('t_i');
% ylabel('t_j');
% zlabel('R_x(t_i, t_j)');
% title('Cross-Correlation Matrix');
% 
% % Optional: Rotate the plot for better visualization
% rotate3d on;
