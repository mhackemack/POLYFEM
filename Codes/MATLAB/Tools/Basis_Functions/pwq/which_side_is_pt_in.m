function out = which_side_is_pt_in(xq,yq,sides)

for k=1:length(sides)
    [in(k),on(k)]=inpolygon(xq,yq,sides{k}.x,sides{k}.y);
end

ind=find(in==1);

howmany_sides = length(ind);

switch howmany_sides
    % inside only 1 side
    case{1}
        out = ind;
        % in between 2 sides
    case{2,3,4}
        if isequal(in,on)
            % we pick the first edge, it doesn't matter which one, the value will be the same
            out=ind(1);
        else
            [in on]
            sides{1}.x
            sides{2}.x
            sides{3}.x
            sides{4}.x
            error('pt %g,%g was found in more than one side %d %d %d %d \n',xq,yq,ind)
        end
    otherwise
        [in on]
        sides{1}.x
        sides{2}.x
        sides{3}.x
        sides{4}.x
        error('pt %g,%g \n',xq,yq)
end

if out<1 || out>4
    error('only 4 sides allowed');
end
