%% Nice Legend ============================================================
% Description: This function calls the legend function, and applies some
% nice settings to make publication quality figures.
%
% Inputs:
%   varargin - Set of labels for the legend
%
% Outputs:
%   hLegend - handle for the legend 
%
% Other Functions Called:
%   legend - the standard MATLAB legend function
%
% Created by: Cory Fraser - Sept 13, 2018
% Latest Edit: Cory Fraser - Sept 13, 2018

%% ========================================================================
function hLegend = NiceLegend(varargin)

labels = [varargin];
for i = 1:length(varargin{:})
    varargin{:}(i) = strcat(varargin{:}(i), ' \hspace{0.25cm} ');
end

hLegend = legend(varargin{:},'orientation','horizontal', 'box','off','color','none',...
        'Position',[0.65 0.9 0.1 0.1], 'Interpreter', 'Latex');

end