%   Example of mesh generation -- square with hole mesh
function [p,t,NIN]=sample_mesh(h)

disp([' ']);
disp(['Star with hole, h=',num2str(h)]);
fdx=@(p) ddiff(dstar(p), dcircle(p,0,0,0.25));
fh= @(p) min(8*sqrt(sum(p.^2,2))-1,2);
box=[-1,-1;1,1];
theta = (0:.1:20)*pi/10;
r = 0.75 + 0.25*sin(5*theta);
x = r.*cos(theta);
y = r.*sin(theta);
theta2 = (0:0.5:20)*pi/10;
x2 = 0.25*cos(theta2);
y2 = 0.25*sin(theta2);
%fix = [x',y'];
fix = [x',y';x2', y2'];

%fix=[-1,-1;-1,1;1,-1;1,1];
[p,t]=distmesh2d(fdx,@huniform,h,box,fix);
%[p,t]=distmesh2d(fdx,@(p) fh(p),h,box,fix);
[p,t]=fixmesh(p,t);
post(p,t,@huniform)

[p,t,NIN]=boundary_reorder(p,t,h,fdx);

figure(1); hold on; plot(p(NIN+1:end,1),p(NIN+1:end,2),'ro'); hold off;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function post(p,t,fh,varargin)

q=simpqual(p,t);
u=uniformity(p,t,fh,varargin{:});
disp(sprintf(' - Min quality %.2f',min(q)))
disp(sprintf(' - Uniformity %.1f%%',100*u))
disp(' ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [h] = fh(p)
%     [theta, rho] = cart2pol(p(:,1),p(:,2));
%     h = min(rho-0.25, .001);
%     h = min(h, rho - 0.75 - 0.25*sin(5*theta));
% end