function varargout = emptyFunction(varargin)
if nargout > 0
    varargout = cell(nargout,1);
    for i=1:nargout
        varargout{i} = 0;
    end
end