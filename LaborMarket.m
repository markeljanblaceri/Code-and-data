%% Labor Market Flows in Albania 
% Estimating Job Finding and Employment Exit Probabilities 
% Author: Markeljan Blacëri
% E-mail: markeljanblaceri@gmail.com

clear;      
clc;        
close all; 

%% Data
filename = 'Albania JFP.xlsx'; 
sheet = 'Data'; 
data = readtable(filename, 'Sheet', sheet);
disp(data(1:5, :)); % Display the first few rows to confirm

% Define the variable
u_r = data.urate % Unemployment rate
u_r_s = data.sr_urate % Short-term unemployment rate
e_l = data.employment % Employment level
lf = data.labor_force % Labor force level 
u = (u_r ./ 100) .* lf % Unemployment level
u_s = (u_r_s ./ 100) .* lf % Short-term unemployment level 

%% Estimating the natural rate of unemployment using the HP filter

% Apply the Hodrick-Prescott filter 
[trend_u, cycle_u] = hpfilter(u_r);

% The trend component (trend_u) represents the estimated natural rate of unemployment
% The cycle component (cycle_u) represents the deviation from the natural rate

% Time settings
start_year = 2012;
num_quarters = length(u_r); % Ensure this matches the length of u_r

% Create a time vector with quarter labels
quarters = cell(1, num_quarters); 
for i = 1:num_quarters
    year = start_year + floor((i - 1) / 4);  
    quarter = mod(i - 1, 4) + 1;             
    quarters{i} = sprintf('%d Q%d', year, quarter); 
end

% Plot the original unemployment series and the estimated natural rate (trend)
figure; 
hold on;

% Plot original unemployment data (u_r)
plot(1:length(u_r), u_r, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Unemployment Rate');

% Plot the estimated natural rate of unemployment (trend_u)
plot(1:length(trend_u), trend_u, '-s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Natural Rate of Unemployment');

hold off; % Release the plot hold

% Customize plot
xlabel('Quarter');
ylabel('Percent');
legend('show');
xticks(1:length(quarters)); 
xticklabels(quarters); 
xtickangle(45);
xlim([1, length(quarters)]);
grid off; 

%% Job finding probability with unemployment rates
% Offset u_r and u_r_s for t+1
u_r_t_plus_1 = u_r(2:end);        % Unemployment rate at t+1
u_r_s_t_plus_1 = u_r_s(2:end);    % Short-term unemployment rate at t+1
u_r_t = u_r(1:end-1);             % Unemployment rate at t

% Compute F_t
F_t = 1 - (u_r_t_plus_1 - u_r_s_t_plus_1) ./ u_r_t;
% Display 
disp(F_t);

% Compute job finding rate f_t = -ln(1-F_t)
f_t = -log(1 - F_t);
% Display 
disp(f_t);

% Plot
start_year = 2012;
end_year = 2023;
num_quarters = length(F_t);

% Create a time vector with quarter labels (e.g., '2012 Q1', '2012 Q2', ..., '2023 Q4')
quarters = cell(1, num_quarters);
for i = 1:num_quarters
    year = start_year + floor((i - 1) / 4);  
    quarter = mod(i - 1, 4) + 1;             
    quarters{i} = sprintf('%d Q%d', year, quarter); 
end

% Plot F_t
figure; 
plot(1:num_quarters, F_t, '-o', 'LineWidth', 2, 'MarkerSize', 6); 
grid off; 
xlabel('Quarter'); 
ylabel('Job finding probability'); 

% Customize x-axis with quarter labels
xticks(1:num_quarters);            
xticklabels(quarters);             
xtickangle(45);                    
xlim([1 num_quarters]);
set(gca, 'Box', 'off');

%% Job finding probability with number of unemployed individuals
% Offset u and u_s for t+1
u_t_plus_1 = u(2:end);        % Unemployment level at t+1
u_s_t_plus_1 = u_s(2:end);    % Short-term unemployment level at t+1
u_t = u(1:end-1);             % Unemployment level at t

% Compute F_t
F_t_new = 1 - (u_t_plus_1 - u_s_t_plus_1) ./ u_t;
% Display 
disp(F_t_new);

% Compute job finding rate f_t = -ln(1-F_t)
f_t_new = -log(1 - F_t_new);
% Display 
disp(f_t_new);

%% Plotting the job finding probabilities
start_year = 2012;
end_year = 2023;
num_quarters = length(F_t);

quarters = {};
for i = 1:num_quarters
    year = start_year + floor((i - 1) / 4);  
    quarter = mod(i - 1, 4) + 1;            
    quarters{i} = sprintf('%d Q%d', year, quarter);
end

% Plot both F_t and F_t_new
figure;
plot(F_t, '-o', 'LineWidth', 2, 'DisplayName', 'F_t^{rate}');
hold on;
plot(F_t_new, '-s', 'LineWidth', 2, 'DisplayName', 'F_t');
hold off;


xlabel('Quarter');
ylabel('Job Finding Probability');
legend('show');
xticks(1:num_quarters);       
xticklabels(quarters);       
xtickangle(45);               
xlim([1, num_quarters]);
grid off;
set(gca, 'Box', 'off');

%% Employment exit rate

n = length(u); % Number of time steps
x_t = zeros(1, n); % Initialize x_t array

% Define a tolerance for convergence
tolerance = 1e-6;

% Iterate over the time series
for t = 1:n-1 
    g = @(x) ((1 - exp(-f_t_new(t) - x)) * x / (f_t_new(t) + x) * lf(t) + exp(-f_t_new(t) - x) * u(t)) - u(t + 1);

    % Solve for x_t using fsolve
    options = optimset('Display', 'off', 'TolFun', tolerance);
    x_t(t) = fsolve(g, max(0, u(t) / lf(t)), options); % Initial guess based on data
    
    % Ensure x_t >= 0
    if x_t(t) < 0
        x_t(t) = 0;
    end
end

% Display results
disp('Computed x_t values:');
disp(x_t);

%% Calculate employment exit probability 
% Calculate X_t from x_t
X_t = 1 - exp(-x_t);

% Display results
disp('Computed X_t values:');
disp(X_t);

%% Employment exit probability - alternative approach
n = length(u); 
X_t_new = zeros(1, n-1); 

% Iterate over the time series
for t = 1:n-1
    X_t_new(t) = (u(t+1) - (1 - F_t_new(t)) * u(t)) / e_l(t);
end

% Display results
disp('Computed X_t_new values:');
disp(X_t_new);

%% Plot
start_year = 2012;
end_year = 2023;

% Number of quarters (assuming num_quarters matches your data length)
num_quarters = length(X_t);

X_t_plot = X_t(1:end-1); 
X_t_new_plot = X_t_new;   

% Create a time vector with quarter labels
quarters = cell(1, num_quarters-1); 
for i = 1:num_quarters-1
    year = start_year + floor((i - 1) / 4);  
    quarter = mod(i - 1, 4) + 1;             
    quarters{i} = sprintf('%d Q%d', year, quarter); 
end

% Plot X_t and X_t_new values
figure; 
hold on;

% Plot X_t
plot(1:num_quarters-1, X_t_plot, '-o', 'LineWidth', 1.5, 'MarkerSize', 6); 

% Plot X_t_new
plot(1:num_quarters-1, X_t_new_plot, '-s', 'LineWidth', 1.5, 'MarkerSize', 6); 

hold off; 

% Customize plot
xlabel('Quarter', 'Interpreter', 'latex'); 
ylabel('Employment exit probability', 'Interpreter', 'latex'); 
xticks(1:num_quarters-1); 
xticklabels(quarters); 
xtickangle(45); 
xlim([1 num_quarters-1]); 
grid off;
legend({'$X_t$', '$\tilde{X}_t$'}, 'Interpreter', 'latex'); 
set(gca, 'Box', 'off'); 