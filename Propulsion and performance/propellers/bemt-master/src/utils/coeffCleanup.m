function dataOut = coeffCleanup(dataIn, rpl, lower, upper)
%COEFFCLEANUP removes negative and exploding values then replaces them
%   dataIn: data to ckean up
%   rpl: optional parameter to choose what the impossible values are
%   replace with (default = 0)
%
%	Usage:
%       dataOut = coeffCleanup(dataIn)
%       dataOut = coeffCleanup(dataIn, rpl, lower, upper)
%       dataOut = coeffCleanup(dataIn, rpl, lower, upper)
%
%

% optional parameters
if nargin == 1
    rpl = 0;
    lower = 0;
    upper = 1;
end

if nargin == 2
    lower = 0;
    upper = 1;
end

% Cleanup
selected_val = (dataIn < lower | dataIn > upper);   % finding impossible values

dataOut = dataIn .* (~selected_val);        % all selected values are set to 0
dataOut = dataOut + selected_val * rpl;     % the replacement is added where needed

% first value if second is replaced
if dataOut(1) > lower && dataOut(2) == rpl
    dataOut(1) = rpl;
end

% isolated values
for i = 2:1:length(dataIn)-1
    if dataOut(i) > lower && dataOut(i+1) == rpl && dataOut(i-1) == rpl
        dataOut(i) = rpl;
    end
end

end

