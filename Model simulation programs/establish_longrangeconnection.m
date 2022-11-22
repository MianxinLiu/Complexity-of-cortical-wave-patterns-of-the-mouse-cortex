linkout_site=[];
receive_site=[];


% 1. region-to-region
% for i1=1:30
%     for j1=1:30
%         linkout_site=[linkout_site;[i1,j1]];
%     end
% end
% for i1=51:80
%     for j1=51:80
%         receive_site=[receive_site;[i1,j1]];
%     end
% end

% 2. point-to-point
linkout_site=[63,38];
receive_site=[38,63];

LRposition.linkout=linkout_site;
LRposition.receive=receive_site;



