% averef() - convert common-reference EEG data to average reference
%
% Usage:
%   >> data = averef(data);
%   >> [data_out W_out S_out meandata] = averef(data,W);
%
% Inputs:
%   data - 2D data matrix (chans,frames*epochs) 
%   W    - ICA weight matrix
%
% Outputs:
%   data_out - Input data converted to average reference.
%   W_out    - ICA weight matrix converted to average reference
%   S_out    - ICA sphere matrix converted to eye()
%   meandata - (1,dataframes) mean removed from each data frame (point)
%
% Note: If 2 args, also converts the weight matrix W to average reference:
%         If ica_act = W*data, then data = inv(W)*ica_act; 
%         If R*data is the average-referenced data, 
%         R*data=(R*inv(W))*ica_act and W_out = inv(R*inv(W));
%         The average-reference ICA maps are the columns of inv(W_out).
%
% Authors: Scott Makeig and Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 1999 
%
% See also: reref()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: averef.m,v $
% Revision 1.9  2002/11/15 01:40:38  arno
% header for web
%
% Revision 1.8  2002/11/12 18:01:21  arno
% header typo
%
% Revision 1.7  2002/09/05 00:30:23  scott
% added meandata output -sm
%
% Revision 1.6  2002/08/21 02:08:19  arno
% nothing
%
% Revision 1.5  2002/08/21 02:03:32  arno
% debugging ica averef
%
% Revision 1.4  2002/08/21 00:21:51  arno
% debugging
%
% Revision 1.3  2002/04/11 18:37:33  scott
% revised help msg
%
% Revision 1.2  2002/04/11 18:02:03  arno
% computing average reference of components
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 12/16/99 Corrected denomiator on the suggestion of Ian Nimmo-Smith, Cambridge UK
% 01-25-02 reformated help & license -ad 

function [data, W, S, meandata] = averef(data, W, S)

if nargin<1
  help averef
  return
end
chans = size(data,1);
if chans < 2 
  help averef
  return
end

% avematrix = eye(chans)-ones(chans)*1/chans;
% data = avematrix*data; % implement as a matrix multiply
% else (faster?)

meandata = sum(data)/chans;
data = data - ones(chans,1)*meandata;

% treat optional ica parameters
if nargin == 2
	winv  = pinv(W);
    size1 = size(winv,1);
	avematrix = eye(size1)-ones(size1)*1/size1;
	W = pinv(avematrix*winv);
end;
if nargin >= 3
	winv = pinv(W*S);
    size1 = size(winv,1);
	avematrix = eye(size1)-ones(size1)*1/size1;
	W = pinv(avematrix*winv);
	S = eye(chans);
end;
