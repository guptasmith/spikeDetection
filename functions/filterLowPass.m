function filteredData = filterLowPass(dataRow,freq,newparams)
% filteredData = filterLowPass(dataRow,freq,newparams). 
% Low-pass filters a time series (one Row) using order-3 Butterworth filter
% cutoff frequency must be provided in the second argument
% Default parameters
% 	params.Fs = 15000;
% Outputs filtered data

% Default parameters
params.Fs = 15000;

% update the parameters 
if nargin > 2
    params = updateParams(params,newparams);
end


Nyquist_frequency=params.Fs/2;
[b,  a]=butter(3,  freq/Nyquist_frequency, 'low');
filteredData=filtfilt(b,  a,  dataRow);

return