function md = adjust_thk_ols(md)
    assert(length(md.geometry.bed) == md.mesh.numberofvertices);
    md.geometry.surface = nan*ones(md.mesh.numberofvertices, 1);
    md.geometry.base = nan*ones(md.mesh.numberofvertices, 1);

    floating_base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness;
    pos = find(floating_base>md.geometry.bed);
    ground_pos = find(floating_base<=md.geometry.bed);
    md.geometry.base(pos) = floating_base(pos);
    md.geometry.base(ground_pos) = md.geometry.bed(ground_pos);
    md.geometry.surface = md.geometry.base + md.geometry.thickness;
    di=md.materials.rho_ice/md.materials.rho_water;
    md.mask.ocean_levelset=md.geometry.thickness+md.geometry.bed/di;
end