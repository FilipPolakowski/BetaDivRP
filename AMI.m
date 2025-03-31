function [taus, ami_values] = AMI(time_series, max_tau, num_bins)
    % Compute Auto-Mutual Information (AMI) for a given time series.
    %
    % Parameters:
    % - time_series: 1D array of time series data
    % - max_tau: Maximum time delay to compute AMI
    % - num_bins: Number of bins for probability estimation
    %
    % Returns:
    % - taus: Array of delay values
    % - ami_values: AMI values for corresponding delays

    taus = 1:max_tau;
    ami_values = zeros(size(taus));

    % Discretize the time series into bins
    [counts, edges] = histcounts(time_series, num_bins);
    
    % Digitize the time series
    digitized_series = discretize(time_series, edges);

    for i = 1:length(taus)
        tau = taus(i);
        x = digitized_series(1:end-tau);   % Original series
        y = digitized_series(tau+1:end);   % Delayed series

        % Compute mutual information
        ami_values(i) = mutual_info(x, y, num_bins);
    end
end

function mi = mutual_info(x, y, num_bins)
    % Compute mutual information between two discrete variables
    joint_hist = accumarray([x(:), y(:)], 1, [num_bins, num_bins]);
    joint_prob = joint_hist / sum(joint_hist(:));

    % Compute marginal probabilities
    p_x = sum(joint_prob, 2);
    p_y = sum(joint_prob, 1);

    % Compute mutual information
    [rows, cols] = find(joint_prob > 0);
    mi = sum(joint_prob(sub2ind(size(joint_prob), rows, cols)) .* ...
        log2(joint_prob(sub2ind(size(joint_prob), rows, cols)) ./ ...
        (p_x(rows) .* p_y(cols)')));
end

