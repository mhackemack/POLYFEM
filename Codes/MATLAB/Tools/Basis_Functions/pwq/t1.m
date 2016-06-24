function val = t1(xq,yq,sides)

side_ID = which_side_is_pt_in(xq,yq,sides);

switch side_ID
    case{4,1}
        val=1-xq-yq;
    case{2,3}
        val=0;
end

