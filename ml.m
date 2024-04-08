% Data generation
num_samples = 1000;
input_data = zeros(num_samples, 3); % [valence, np_radius, koff]
output_data = zeros(num_samples, 1); % binding frequency

for i = 1:num_samples
    vh = randi([1, 100]);
    np_radius = randi([50, 500]);
    koff = 0.01 + (10-0.01).*rand(1,1); % random koff between 0.01 and 10
    
    np_params = [vh, np_radius, 1];
    kinetic_params = [0.01, 0.1, koff];
    tcr_params = [1000, 15, 50, 20, 300]; % tcr_params_cluster
    
    bind_freq = NP_surface_binding(np_params, kinetic_params, tcr_params, 0);
    
    input_data(i, :) = [vh, np_radius, koff];
    output_data(i) = bind_freq;
end
% Train a neural network regression model

% Split data into training and test sets
train_ratio = 0.8;
[trainInd, ~, testInd] = dividerand(num_samples, train_ratio, 0, 1 - train_ratio);

x_train = input_data(trainInd, :)';
y_train = output_data(trainInd)';

x_test = input_data(testInd, :)';
y_test = output_data(testInd)';

% Define neural network
net = feedforwardnet([10, 10]); % Two hidden layers with 10 neurons each
net.trainFcn = 'trainlm'; % Levenberg-Marquardt backpropagation.
net.divideFcn = 'divideind';
net.divideParam.trainInd = trainInd;
net.divideParam.testInd = testInd;

% Train the network
net = train(net, x_train, y_train);

% Evaluate on test set
predictions = net(x_test);
mse_error = mse(y_test - predictions);

fprintf('Mean Squared Error on Test Set: %f\n', mse_error);
% Assuming you already have trained the neural network as 'net'
% and your test data as x_test and y_test from the previous code.

% Compute predictions on original test data
original_predictions = net(x_test);
original_mse = mse(y_test - original_predictions);

% Number of features
num_features = size(x_test, 1);

% Initialize array to store feature importances
feature_importances = zeros(num_features, 1);

for i = 1:num_features
    % Copy the original test data
    permuted_x_test = x_test;
    
    % Shuffle the i-th feature
    permuted_x_test(i, :) = permuted_x_test(i, randperm(size(permuted_x_test, 2)));
    
    % Compute predictions with the permuted data
    permuted_predictions = net(permuted_x_test);
    permuted_mse = mse(y_test - permuted_predictions);
    
    % Feature importance is the drop in performance
    feature_importances(i) = permuted_mse - original_mse;
end

% Display feature importances
fprintf('Feature importance for valence: %f\n', feature_importances(1));
fprintf('Feature importance for np_radius: %f\n', feature_importances(2));
fprintf('Feature importance for koff: %f\n', feature_importances(3));
