function val = tc(xq,yq,sides)

side_ID = which_side_is_pt_in(xq,yq,sides);

switch side_ID
    case{1}
        val=2*yq;
    case{2}
        val=2*(1-xq);
    case{3}
        val=2*(1-yq);
    case{4}
        val=2*xq;
end

