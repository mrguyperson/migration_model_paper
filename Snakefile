rule targets:
    input:
        "data/ds2890.zip",
        "data/ds2890.gdb",
        "data/cover.shp"
        
rule download_archive:
    input:
        script = "code/get_shape_files.bash"
    output:
        "data/ds2890.zip"
    shell:
        """
        {input.script}
        """

rule unzip_archive:
    input:
        script = "unzip_shape_files.bash"
    output:
        "data/ds2890.gdb"
    shell:
        """
        {input.script}
        """
rule convert_to_shapefile:
    input:
        script = "simplify_shape_file.R"
    output:
        "data/cover.shp"
    shell:
        """
        {input.script}
        """
