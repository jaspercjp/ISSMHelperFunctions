
% calc_area_from_elt(md, elts) uses the shoelace formula to calculate the total area
% encased by all the triangles given in the elts array.
function elt_area=calc_area_from_elt(md, elts)
    elt_area=0;
    for i=1:length(elts)
        x=md.mesh.x(elts{i});
        y=md.mesh.y(elts{i});
        elt_area=elt_area + abs(x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2)))/2;
    end
end