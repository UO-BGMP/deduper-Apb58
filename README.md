# PCR Duplicate Read Remover

#### Adrian Bubie
#### 1/1/2018

### Summary:
This program removes putative duplicate reads from a given SAM file that are the result of PCR amplification during sequencing library prep. A new, "deduped" version of the SAM file is produced in the same directory as the input SAM file. Several arguments can be specified when running the script (see below). 

There is no installation required to use this script, but the machine running the program must have `python3.4+`.

### Running the Script:
Download `PCR_Deduper.py` and edit the path to point to your local version of python3 (see `#!/usr/bin/python3` at the top of the file). You will also need to make the script executable, which can be done from the command line with:

```
$ chmod 777 PCR_Deduper.py
```

To run from the command line (Linux/UNIX/MacOS), navigate to the directory where the script is located and use:

```
$ ./PCR_Deduper.py -f {path/to/Sam_file.sam}
```

This will run the script with the default options on the given SAM file, and write the resulting duplicate-less SAM to the same directory.

For additional running options, or help, see `./PCR_Deduper.py -h`.
