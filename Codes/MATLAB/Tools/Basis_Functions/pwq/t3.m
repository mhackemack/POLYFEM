function val = t3(xq,yq,sides)

side_ID = which_side_is_pt_in(xq,yq,sides);

switch side_ID
    case{2,3}
        val=xq+yq-1;
    case{4,1}
        val=0;
end

