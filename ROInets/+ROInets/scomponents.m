function [sci, sizes] = scomponents(A)
% SCOMPONENTS Compute the strongly connected components of a graph
%
% ci=scomponents(A) returns an index for the component number of every 
% vertex in the graph A.  The total number of components is max(ci).
% If the input is undirected, then this algorithm outputs just the 
% connected components.  Otherwise, it output the strongly connected
% components.
%
% The implementation is from Tarjan's 1972 paper: Depth-first search and 
% linear graph algorithms. In SIAM's Journal of Computing, 1972, 1, 
% pp.146-160.
%
% See also DMPERM
%
% Example:
%   load_gaimc_graph('cores_example'); % the graph A has three components
%   ci = scomponents(A)
%   ncomp = max(ci)               % should be 3
%   R = sparse(1:size(A,1),ci,1,size(A,1),ncomp); % create a restriction matrix
%   CG = R'*A*R;                  % create the graph with each component 
%                                 % collapsed into a single node.

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2008-04-11: Initial coding

if isstruct(A), rp=A.rp; ci=A.ci; %ofs=A.offset;
else [rp, ci]=sparse_to_csr(A); 
end

n=length(rp)-1; sci=zeros(n,1); cn=1;
root=zeros(n,1); dt=zeros(n,1); t=0;
cs=zeros(n,1); css=0; % component stack
rs=zeros(2*n,1); rss=0; % recursion stack holds two nums (v,ri)
% start dfs at 1
for sv=1:n
    v=sv; if root(v)>0, continue; end
    rss=rss+1; rs(2*rss-1)=v; rs(2*rss)=rp(v); % add v to the stack
    root(v)=v; sci(v)=-1; dt(v)=t; t=t+1;
    css=css+1; cs(css)=v; % add w to component stack
    while rss>0
        v=rs(2*rss-1); ri=rs(2*rss); rss=rss-1; % pop v from the stack
        while ri<rp(v+1)
            w=ci(ri); ri=ri+1;
            if root(w)==0
                root(w)=w; sci(w)=-1;  dt(w)=t; t=t+1;
                css=css+1; cs(css)=w; % add w to component stack
                rss=rss+1; rs(2*rss-1)=v; rs(2*rss)=ri; % add v to the stack
                v=w; ri=rp(w); continue; 
            end
        end
        for ri=rp(v):rp(v+1)-1
            w=ci(ri); 
            if sci(w)==-1
                if dt(root(v))>dt(root(w)), root(v)=root(w); end
            end
        end
        if root(v)==v
            while css>0
                w=cs(css); css=css-1; sci(w)=cn;
                if w==v, break; end
            end
            cn=cn+1;
        end
    end
end

if nargout>1
    sizes=accumarray(sci,1,[max(sci) 1]);
end
end%scomponents

function [rp, ci, ai, ncol]=sparse_to_csr(A,varargin)
% SPARSE_TO_CSR Convert a sparse matrix into compressed row storage arrays
% 
% [rp ci ai] = sparse_to_csr(A) returns the row pointer (rp), column index
% (ci) and value index (ai) arrays of a compressed sparse representation of
% the matrix A.
%
% [rp ci ai] = sparse_to_csr(i,j,v,n) returns a csr representation of the
% index sets i,j,v with n rows.
%
% Example:
%   A=sparse(6,6); A(1,1)=5; A(1,5)=2; A(2,3)=-1; A(4,1)=1; A(5,6)=1; 
%   [rp ci ai]=sparse_to_csr(A)
%
% See also CSR_TO_SPARSE, SPARSE

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2008-04-07: Initial version
% 2008-04-24: Added triple array input
% 2009-05-01: Added ncol output
% 2009-05-15: Fixed triplet input

narginchk(1, 5)
retc = nargout>1; reta = nargout>2;

if nargin>1
    if nargin>4, ncol = varargin{4}; end
    nzi = A; nzj = varargin{1};
    if reta && length(varargin) > 2, nzv = varargin{2}; end    
    if nargin<4, n=max(nzi); else n=varargin{3}; end
    nz = length(A);
    if length(nzi) ~= length(nzj), error('gaimc:invalidInput',...
            'length of nzi (%i) not equal to length of nzj (%i)', nz, ...
            length(nzj)); 
    end
    if reta && length(varargin) < 3, error('gaimc:invalidInput',...
            'no value array passed for triplet input, see usage'); 
    end
    if ~isscalar(n), error('gaimc:invalidInput',...
            ['the 4th input to sparse_to_csr with triple input was not ' ...
             'a scalar']); 
    end
    if nargin < 5, ncol = max(nzj); 
    elseif ~isscalar(ncol), error('gaimc:invalidInput',...
            ['the 5th input to sparse_to_csr with triple input was not ' ...
             'a scalar']); 
    end
else
    n = size(A,1); nz = nnz(A); ncol = size(A,2);
    retc = nargout>1; reta = nargout>2;
    if reta,     [nzi, nzj, nzv] = find(A); 
    else         [nzi, nzj]      = find(A);
    end
end
if retc, ci = zeros(nz,1); end
if reta, ai = zeros(nz,1); end
rp = zeros(n+1,1);
for i=1:nz
    rp(nzi(i)+1)=rp(nzi(i)+1)+1;
end
rp=cumsum(rp);
if ~retc && ~reta, rp=rp+1; return; end
for i=1:nz
    if reta, ai(rp(nzi(i))+1)=nzv(i); end
    ci(rp(nzi(i))+1)=nzj(i);
    rp(nzi(i))=rp(nzi(i))+1;
end
for i=n:-1:1
    rp(i+1)=rp(i);
end
rp(1)=0;
rp=rp+1;
end%sparse_to_csr
% [EOF]