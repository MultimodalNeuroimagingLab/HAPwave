function events_table = bids_clipEvents(events_table,key,val)

%   >> mefObj.filterEvents(Name, Value, ...);
%       Name =                  char array. Column name in mefObj.evts to indicate which column to look for matches. E.g. 'electrical_stimulation_current'
%       Value =                 char array or cell array of chars to match in mefObj.evts.<Name>. E.g. '6.0 mA' or {'4.0 mA', '6.0 mA'}
% 
% usage example
% events_table_clipped = bids_clipEvents(events_table,'electrical_stimulation_current', {'4.0 mA', '6.0 mA'}); % keep only events with stim current == 4.0 or 6.0 mA

valStr = join(val, ', ');
if ~sum(ismember(events_table.(key), val))
    error('Cannot find any matches for %s in obj.%s', valStr{1}, key);
end
events_table(~ismember(events_table.(key), val), :) = [];

