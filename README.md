# PCR Duplicate Read Remover

#### Adrian Bubie
#### 1/10/2018

### Summary:
This program removes putative duplicate reads from a given SAM file that are the result of PCR amplification during sequencing library prep. A new, "deduped" version of the SAM file is produced in the same directory as the input SAM file. Several arguments can be specified when running the script (see below). 

There is no installation required to use this script, but the machine running the program must have `python3.4+`.

**NOTE:** There are currently 2 versions of this script; `PCR_Deduper_v1.py` is less memory intensive upfront, but is much slower for large (>500000) readsets. `PCR_Deduper_v2.py` usesa different approach than described in the documentation, and requires memory enough to read/write in your entire read set in duplicate, which can be a limitation for VERY (>2GB) SAM files. However, it is *much* faster than v1, so if the memory is not an issue, v2 should be preferred.

<br>

### Running the Script:
Download `PCR_Deduper_v2.py` and edit the path to point to your local version of python3 (see `#!/usr/bin/python3` at the top of the file). You will also need to make the script executable, which can be done from the command line with:

```
$ chmod 777 PCR_Deduper_v2.py
```

To run from the command line (Linux/UNIX/MacOS), navigate to the directory where the script is located and use:

```
$ ./PCR_Deduper_v2.py -f {path/to/Sam_file.sam}
```

This will run the script with the default options on the given SAM file, and write the resulting duplicate-less SAM to the same directory.

For additional running options, or help, see `./PCR_Deduper_v2.py -h`.
