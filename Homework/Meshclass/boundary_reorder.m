function [pnew,tnew,ke]=boundary_reorder(p,t,h0,fdx)

N=length(p);

ks=1;
ke=N;
mb=0;

pnew = p;
tnew = t;

geps=0.001*h0;

for i=1:N
    if fdx(p(i,1:2))<-geps
        pnew(ks,:) = p(i,:);
        ks = ks + 1;
    else
        pnew(ke,:)=p(i,:);
        ke = ke - 1;
    end
end

tnew=delaunayn(pnew);                            % re-triangulate
pmid=(pnew(tnew(:,1),:)+pnew(tnew(:,2),:)+pnew(tnew(:,3),:))/3;    % Compute centroids
tnew=tnew(feval(fdx,pmid)<-geps,:);         % Keep interior triangles
%[pnew,tnew]=fixmesh(pnew,tnew);