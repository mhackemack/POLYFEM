function val = t2(xq,yq,sides)

side_ID = which_side_is_pt_in(xq,yq,sides);

switch side_ID
    case{1,2}
        val=xq-yq;
    case{3,4}
        val=0;
end

