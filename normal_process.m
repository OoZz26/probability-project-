% Set the parameters
num_realizations = 1000;
t = 0:0.01:2; % Time vector from 0 to 2 with step 0.1
numberofelements=0:1:999;

% Initialize the matrix to store realizations
X = zeros(num_realizations, length(t));

% Generate realizations
A = normrnd(-5,5,1,num_realizations); % Draw A from a normal distribution with mean 55
b=A';
disp(b(1,:));
for i = 1:num_realizations

    W_t = b(i,:) * cos(4 * t);
    X(i, :) = W_t;
end
% mean(X(:,1))
% Print or use the matrix as needed
%plot(t,X(1, :))
% acf1=zeros(num_realizations,401);
% for m =1:num_realizations
%    acf1(m,:) = xcorr(X(m,:)); % 'biased' option normalizes by N, 'unbiased' normalizes by N-k
% end
% % % Calculate autocorrelation function
% % for m =1:num_realizations
% % 
% %     for k = 0:length(t)-1
% %         for n = 1:length(t)-k
% %             acf1(m,k+1) = acf1(m,k+1) + X(m,n) * X(m,n+k);
% %         end
% %         acf1(m,k+1) = acf1(m,k+1) / length(t); % Normalize by the total number of samples
% %     end
% % end
% % % Normalize the ACF at zero lag to 1
% % acf1 = acf1 ./ acf1(1,1);
% 
% % Plot the autocorrelation function
% % figure;
% surf(numberofelements,numberofelements,acf1)
