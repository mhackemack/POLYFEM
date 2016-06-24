function val = t4(xq,yq,sides)

side_ID = which_side_is_pt_in(xq,yq,sides);

switch side_ID
    case{3,4}
        val=yq-xq;
    case{1,2}
        val=0;
end

