function md = delete_prev_transientsol(md)
    allFields = fieldnames(md.results);
    fieldToKeep = 'TransientSolution';
    fieldsToDelete = allFields(~strcmp(allFields, fieldToKeep));
    md.results = rmfield(md.results, fieldsToDelete);
end