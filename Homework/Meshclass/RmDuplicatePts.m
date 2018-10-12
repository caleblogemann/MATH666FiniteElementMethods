function [p]=RmDuplicatePts(p)

snap=max(max(p,[],1)-min(p,[],1),[],2)*1024*eps;
[foo,ix,jx]=unique(round(p/snap)*snap,'rows');
p=p(ix,:);