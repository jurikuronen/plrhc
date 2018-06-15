function [good tests] = test_data(data, N)
% Tests if the data is in the right format.
% Juri Kuronen (2018)
%
% Input:
% - data: discrete Nxd matrix
% - N: number of data samples (optional)
%
% Output:
% - good: positive value indicates data is in the right format
% - tests: individual test results

    % Outcomes should be {1, 2, ...}
    u = unique(data);
    test1 = 0;
    for i = 1:length(u)
        test1 = (test1 || (u(i) == i));
    end

    % Data matrix should be in 'double' format
    test2 = isa(data, 'double');

    % No infs or nans
    test3 = sum(isnan(data(:))) == 0;
    test4 = sum(isinf(data(:))) == 0;

    tests = [test1 test2 test3 test4];

    % If user supplied intended number of data samples,
    % check that the samples are the rows of the matrix.
    if nargin == 2
        test5 = (N == size(data, 1));
        tests(end + 1) = test5;
    end

    good = all(tests);

end
