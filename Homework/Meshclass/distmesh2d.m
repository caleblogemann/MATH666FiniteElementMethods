function [p,t]=distmesh2d(fd,fh,h0,bbox,pfix,varargin)
%DISTMESH2D 2-D Mesh Generator using Distance Functions.
%   [P,T]=DISTMESH2D(FD,FH,H0,BBOX,PFIX,FPARAMS)
%
%      P:         Node positions (Nx2)
%      T:         Triangle indices (NTx3)
%      FD:        Distance function d(x,y)
%      FH:        Scaled edge length function h(x,y)
%      H0:        Initial edge length
%      BBOX:      Bounding box [xmin,ymin; xmax,ymax]
%      PFIX:      Fixed node positions (NFIXx2)
%      FPARAMS:   Additional parameters passed to FD and FH
%
%   Example: (Uniform Mesh on Unit Circle)
%      fd=inline('sqrt(sum(p.^2,2))-1','p');
%      [p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[]);
%
%   Example: (Rectangle with circular hole, refined at circle boundary)
%      fd=inline('ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5))','p');
%      fh=inline('min(4*sqrt(sum(p.^2,2))-1,2)','p');
%      [p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
%
%   See also: MESHDEMO2D, DISTMESHND, DELAUNAYN, TRIMESH.

%   Copyright (C) 2004-2005 Per-Olof Persson. See COPYRIGHT.TXT for details.

dptol=.001; ttol=.1; Fscale=1.2; deltat=.2; geps=.001*h0; deps=sqrt(eps)*h0;

% 1. Create initial distribution in bounding box (equilateral triangles)
[x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
x(2:2:end,:)=x(2:2:end,:)+h0/2;                      % Shift even rows
p=[x(:),y(:)];                                       % List of node
                                                     % coordinates
p=RmDuplicatePts(p);

% 2. Remove points outside the region, apply the rejection method
p=p(feval(fd,p,varargin{:})<geps,:);                 % Keep only d<0 points
r0=1./feval(fh,p,varargin{:}).^2;                    % Probability to keep point
p=[pfix; p(rand(size(p,1),1)<r0./max(r0),:)];        % Rejection
                                                     % method
% Remove duplicate points, but be careful with fixed points
snap=max(max(p,[],1)-min(p,[],1),[],2)*1024*eps;
[foo,ix,jx]=unique(round(p/snap)*snap,'rows');
p=p(ix,:);
numfixpts = size(pfix,1);
Ifix = jx(1:numfixpts);

N=size(p,1);                                         % Number of points N

pold=inf;                                            % For first iteration
while 1
  % 3. Retriangulation by the Delaunay algorithm
  if max(sqrt(sum((p-pold).^2,2))/h0)>ttol           % Any large movement?
    pold=p;                                          % Save current positions
    t=delaunayn(p);                                  % List of triangles
    pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids
    t=t(feval(fd,pmid,varargin{:})<-geps,:);         % Keep interior triangles
    % 4. Describe each bar by a unique pair of nodes
    bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
    bars=unique(sort(bars,2),'rows');                % Bars as node pairs
    % 5. Graphical output of the current mesh
    figure(1);
    trimesh(t,p(:,1),p(:,2),zeros(N,1));
    set(gcf,'renderer','painters');
    caxis([0 1]);
    d1 = 0.01*(bbox(2,1) - bbox(1,1));
    d2 = 0.01*(bbox(2,2) - bbox(1,2));
    axis([bbox(1,1)-d1 bbox(2,1)+d1 bbox(1,2)-d2 bbox(2,2)+d2]);
    colormap(gray);
    t1 = title('Simple Mesh Generator in MATLAB');
    set(t1,'fontsize',16);
    view(2),axis equal,axis off,drawnow
  end

  % 6. Move mesh points based on bar lengths L and forces F
  barvec=p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
  L=sqrt(sum(barvec.^2,2));                          % L = Bar lengths
  hbars=feval(fh,(p(bars(:,1),:)+p(bars(:,2),:))/2,varargin{:});
  L0=hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));     % L0 = Desired lengths
  F=max(L0-L,0);                                     % Bar forces (scalars)
  Fvec=F./L*[1,1].*barvec;                           % Bar forces (x,y components)
  Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
  Ftot(Ifix,:)=0;                                    % Force = 0 at fixed points
  p=p+deltat*Ftot;                                   % Update node positions

  % 7. Bring outside points back to the boundary
  d=feval(fd,p,varargin{:}); ix=d>0;                 % Find points
                                                     % outside
                                                     % (d>0)
  mstop=0;
  numiter=0;
  while (mstop==0)
      numiter=numiter+1;
      
      dgradx=(feval(fd,[p(ix,1)+deps,p(ix,2)],varargin{:})-d(ix))/deps; % Numerical
      dgrady=(feval(fd,[p(ix,1),p(ix,2)+deps],varargin{:})-d(ix))/deps; % gradient
      dmag = sqrt(dgradx.^2+dgrady.^2);
      dgradx=dgradx./dmag;
      dgrady=dgrady./dmag;
      
      BndProj = [d(ix).*dgradx,d(ix).*dgrady];
      NormBndProj = max(sqrt(BndProj(:,1).^2+BndProj(:,2).^2));
      
      p(ix,:)=p(ix,:)-BndProj;     % Project back to boundary
            
      if (NormBndProj<=1.0e-8)
          mstop = 1;
      else
          d=feval(fd,p,varargin{:});
      end
      if (numiter>1000)
          error('Boundary projection did not terminate');
      end
  end

  % 8. Termination criterion: All interior nodes move less than dptol (scaled)
  d=feval(fd,p,varargin{:});
  if max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0)<dptol, break; end
end
