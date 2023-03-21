rule download_archive:
    input:
        script = "code/get_shape_files.bash"
    output:
        "data/ds2890.zip"
    shell:
        """
        {input.script}
        """