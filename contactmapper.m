function [M,x,y]=contactmapper(s,A,B,varargin)
%%
% 
% CONTACTMAPPER(S,RES)
% CONTACTMAPPER makes all atom contact maps from pdb structs
% between two specified chains
%
% S is s struct output from PDBREAD or filepath to a pdb
% A is the name (string) of the first chain
% B is the name (string) of the second chain
% 
% EXAMPLE CODE:
% 
% s = pdbread('example.pdb')
% [M,x,y] = contactmapper(s,'chain1','chain2)
% imagesc(M)
%
%%
if isempty(regexp(class(s),'struct'))
	s=pdbread(s);
end

if (nargin>3)
	sel=varargin{1};
	ndx=importdata(varargin{2},'[');
	sels=find(~cellfun('isempty',regexp(ndx,'[')));
	if (sel<length(sels))
		atm=str2num(strjoin(ndx(sels(sel)+1:sels(sel+1)-1)'));
	else
		atm=str2num(strjoin(ndx(sels(sel)+1:end)'));
	end
else
	atm=1:length(s.Model.Atom);
end

a=zeros(1,length(atm));
for i=1:length(a)
    a(i)=s.Model.Atom(atm(i)).chainID;
end
a=char(a);


x=atm(regexp(a,A));
y=atm(regexp(a,B));
M=zeros(length(x),length(y));
m=1;
for i=x
	n=1;
    for j=y
        d = [s.Model.Atom(i).X s.Model.Atom(i).Y s.Model.Atom(i).Z] - [s.Model.Atom(j).X s.Model.Atom(j).Y s.Model.Atom(j).Z];
        M(m,n)=sqrt(sum(d.^2));
        n=n+1;
    end
    m=m+1;
end
end