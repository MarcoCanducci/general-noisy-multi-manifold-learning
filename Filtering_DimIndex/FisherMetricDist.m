function Dist = FisherMetricDist(XI,XJ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% FisherMetricDist computes the geodesic distance on simplex between any  %
% two set of points XI and XJ.                                            %
%                                                                         %
%  Copyright (C) 2020  Marco Canducci                                     %
%  Email: marco.canducci91@gmail.com                                      %
%                                                                         %
%  This program is free software: you can redistribute it and/or modify   %
%  it under the terms of the GNU Affero General Public License as         %
%  published by the Free Software Foundation, either version 3 of the     %
%  License, or (at your option) any later version.                        %
%                                                                         %
%  This program is distributed in the hope that it will be useful,        %
%  but WITHOUT ANY WARRANTY; without even the implied warranty of         %
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          %
%  GNU Affero General Public License for more details.                    %
%                                                                         %
% You should have received a copy of the GNU Affero General Public License%
% along with this program.  If not, see <https://www.gnu.org/licenses/>.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Description
%   
%   Input parameters:
%       -   XI: Matrix of objective observations w.r.t. which compute the
%           distances.
%       -   XJ: Matrix of the query obs.
%
%   Output parameters:
%       -   Distance matrix containing the distance to each objective obs
%           in XI of the query obs in XJ.
%
%   The geodesic distances on the d-simplex are computed as the shortest
%   curves onto the positive portion of the d-sphere of radius 2:
%
%   d(x1,x2) = 2*arccos(sum(sqrt(x1*x1')))


XJ = sqrt(XJ);
XI = sqrt(XI);

Dist = acos(XJ*XI');
