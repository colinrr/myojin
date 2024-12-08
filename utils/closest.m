function [oval,ix] = closest(ival,V,N,opt)
% [oval,ix] = closest(ival,vec,opt);
% Find the value and index of closest value in a vector. Vector need not be
% monotonic. Will return N closest, where N by default is 1.
%
% INPUT     ival = value(s) to search for - can be an array of any dimension
%           V    = vector to search for closest values
%           N    = Number of closest values to find.
%
% v NOT IMPLEMENTED v
%           opt  = search option:
%                   'next'=find closest value higher than ival
%                   'prev'=find closest value lower than ival
%                   'abs'= find absolute closest value (default)
%
% OUTPUT    oval = matrix, size [length(ival) x N] of values in vector V 
%                   closest to ival, where each column gives the next
%                   closest value.
%             ix = matrix, size [length(ival) x N] of indices of oval in V
%
%  NOTE: When using 'prev' or 'next' options, it is possible to get
%  redundant indices where a value occurs within N samples of the min/max of vector V.
%
% C Rowell
% Eventually I should set this up for V to be an arbitrary array...

if nargin<3
    N = 1;
end
if nargin<4
    opt='abs';
end

assert(isvector(V),'This function searches vectors. Use closest2d to search a matrix.')
assert(length(V)>=N,'V must have at least N values.')
if isrow(V)
    V = V';
end
if isrow(ival)
    ival = ival';
end

if numel(ival)>1
    ivalSz = size(ival);
    ivalSz = ivalSz(ivalSz>1);
else
    ivalSz = 1;
end
ival = ival(:);

ix   = squeeze(zeros([length(ival) N]));
oval = ix;

[~,VmaxI] = max(V);
[~,VminI] = min(V);

switch opt
    case 'abs'
        
        for ii=1:length(ival)
            [Vmat,I] = sort(abs(ival(ii) - V),'ascend');
            ix(ii,:)   = I(1:N);
        end
        
    case 'next'
        
        for ii=1:length(ival)
            [Vmat,I] = sort(V-ival(ii),'ascend');
            vi = find(Vmat>=0,N,'first');
            ix(ii,1:length(vi))   = I(vi);
        end
        if any(ix(:)==0)
            warning('Index matrix contains empty values - they will be replaced with the index for max(V)')
            ix(ix==0) =  VmaxI;
        end        
        
    case 'prev'
        for ii=1:length(ival)
            [Vmat,I] = sort(V-ival(ii),'ascend');
            vi = flipud(find(Vmat<=0,N,'last'));
            ix(ii,1:length(vi))   = I(vi);
        end 
        if any(ix(:)==0)
            warning('Index matrix contains empty values - they will be replaced with the index for min(V)')
            ix(ix==0) = VminI;
        end        

end

oval = V(ix);


ix   = reshape(ix,[ivalSz N]);
oval = reshape(oval,[ivalSz N]);

    

end