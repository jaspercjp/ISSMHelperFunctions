
% extract_floating_elements(md, ocean_mask, ice_mask) returns a list of all the elements 
% as a cell array of 3x1 vertex arrays.
function float_elt=extract_floating_elements(md, ocean_mask, ice_mask)
    float_vert=find(ocean_mask<0 & ice_mask<0);
    float_elt={};
    count=0;
    for i=1:md.mesh.numberofelements
        elt=md.mesh.elements(i,:);
        if any(float_vert==elt(1)) && any(float_vert==elt(2)) && any(float_vert==elt(3))
            count=count+1;
            float_elt{count}=elt;
        end
    end
end
