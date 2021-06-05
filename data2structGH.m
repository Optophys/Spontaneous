function [data]=data2structGH(cfg)
% This function takes raw LFP data together with SF (obtained from
% and returns a struct usable by field trip
% syntax: [data]=data2structGH(cfg)
%
% Mandatory input:
% ================
% cfg.raw_data      = matrix of channels x timestamps
%
% Optional input:
% ===============
% cfg.SF            = double, the sampling frequency in Hz. Default= 976.5625;
%
% cfg.ch_names      = cell-array of strings with user defined names of the
%                     output channel names ('labels'). 
%                     Default = {'CH001','CH002'...}
% 
% cfg.demean        = 'yes' or 'no', whether to set the mean of the data to
%                     0 (in ft_preprocessing). Default = 'yes'.
% 
% cfg.detrend       = 'yes' or 'no. Default = 'no'.
% 
% cfg.resample      = bool, 0 means no resampling. Default = 0 (don't
%                     resample).
%
% cfg.resamplefs    = double, to which frequency resampling should be done.
%                     Default: the closest kHz to the original sf.
%
% Last updated 01/10/2017, by Golan Karvat

try TDT = cfg.raw_data;
catch fprintf ('Problems loading the raw data\n'); return; end
if isfield(cfg,'SF'); SF = cfg.SF; else SF = 976.5625; end %the default is 976.5625
if ~isfield(cfg,'resample'); cfg.resample = 0; end %the default is not to resample. 
if ~isfield(cfg,'resamplefs'); cfg.resamplefs = round(SF/1000)*1000; end %the default resampling is to the nearest kHz.
if isfield(cfg,'ch_names')
    data.label = cfg.ch_names;
else
    for i=1:size(TDT,1)
        if i<10
            data.label{i,1}=['CH00' num2str(i)];
        else
            if i<100
                data.label{i,1}=['CH0' num2str(i)];
            else
                data.label{i,1}=['CH' num2str(i)];
            end
        end
    end
end
data.fsample=SF;
data.time{1,1}=[1:size(TDT,2)]/SF;
data.trial{1,1}=TDT;
data.nTrials=1;
data.nSamples=size(TDT,2);
data.hdr.FirstTimeStamp=0;

if cfg.resample
    %Resampling      
    if ~isfield(cfg,'detrend'); cfg.detrend    = 'no'; end
    if ~isfield(cfg,'demean'); cfg.demean = 'yes'; end
    [data] = ft_resampledata(cfg, data);
%     data.fsample = cfg.resamplefs;
    data.nSamples = length(data.trial{1});
end
% ensure the data is valid in FieldTrip:
data = ft_preprocessing([],data);