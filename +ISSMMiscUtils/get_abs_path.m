function abs_path = get_abs_path(rel_path)
    javaFile = java.io.File(rel_path);
    abs_path = char(javaFile.getCanonicalPath());
end