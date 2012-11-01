function log = isafq(afq)
% Check if a variable is an afq structure
%
% log = isafq(afq)
%
% Returns a 0 (false) if the input is not an afq structure and a 1 (true)
% if it is

log = isstruct(afq) && isfield(afq,'type') && strcmp(afq.type,'afq version 1');