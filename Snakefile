rule targets:
    input:
        "data/ds2890.zip",
        "data/ds2890.gdb",
        "data/cover.shp",
        "data/cover.prj",
        "data/cover.shx",
        "data/cover.dbf"
        
rule download_archive:
    input:
        script = "code/get_shape_files.bash"
    output:
        "data/ds2890.zip"
    conda:
        "envs/wget.yml"
    container:
        None
    shell:
        """
        {input.script}
        """

rule unzip_archive:
    input:
        script = "code/unzip_shape_files.bash",
        archive = "data/ds2890.zip"
    output:
        directory("data/ds2890.gdb")
    conda:
        "envs/unzip.yml"
    container:
        None
    shell:
        """
        {input.script}
        """
rule convert_to_shapefile:
    input:
        script = "code/simplify_shape_file.R"
    output:
        "data/cover.shp",
        "data/cover.prj",
        "data/cover.shx",
        "data/cover.dbf"
    container: 
        "docker://mrguyperson/fedora36r"
    shell:
        """
        {input.script}
        """
