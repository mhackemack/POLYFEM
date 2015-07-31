function varargout = fixup_node_ordering(~, cnodes, fnodes)
n = nargout;
ind = zeros(1,length(fnodes));
ind2 = zeros(1,length(fnodes));
for i=1:length(fnodes)
    for j=1:length(cnodes)
        if fnodes(i) == cnodes(j)
            ind(i) = cnodes(j);
            ind2(i) = j;
            break
        end
    end
end
if n==1
    varargout{1} = ind;
elseif n==2
    varargout{1} = ind;
    varargout{2} = ind2;
end