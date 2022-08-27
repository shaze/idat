# idat


Simple code to manage idat files


# plotcluster.py


### Requirements

This needs Python3 with pandas, numpy and matplotlib.

### Basic use


This program plots a SNP's intensity file(s). The desired input is  a file which as the Red/Gren intensities. An input file is a text file which has the calls from a number of individuals. Each line is the call of one individual and contains the red/green intensities separated by white space.

The simplest run would be

```plotcluster.py   snp.csv  --output show```

or

```plotcluster.py  snp.csv --output snp.jpg```

This will produce a plot of the given SNP. In the first case, the picture will be shown interactively. You must be running your program locally or, if remotely, using X11 or some other graphical interface.  To send the output to file (either because you want to or have to) give a file name. Note the suffix is compulsory: in the example we created a JPEG file, but you can specify `.png` or `.pdf` if you want.

By default the program does a polar tranformation of the intensities. But note that no normalisation is done so this should only be used for a basic eyeballing test of the clusters.  If you prefer not to do a polar transform you can also specify the `--cartesian` option.

### Basic options

* `--size` : by default the dots in the plots are 5pt in size. You can change this using `--size`
* `--colours` :  to change the default colours you can provide a comma-separated (no spaces) list of colours. For example `--colours pink,cyan,black`
* `--output` :  this is compulsory. Use `--output show` for interactive use or specify the name of the output file with graphics extension (e.g. `--output output.pdf`)


### Plotting multiple intensity files

A common use will be to plot the _same_ SNP calls from different batches. This will allow you to see if there are any marked differences between the two calls. For example, suppose I have two batches that have been genotyped separately and I have extracted out the intensities separately from the two sets of IDAT files. I can plot the the two files like this

```
python3 plotcluster.py batch01/kgp10065099.csv batch02/kgp10065099.csv --output kgp10065099.jpg
```

The calls for batch01 will appear in red, the calls for batch02 in blue. You can then do an eyeball test for obvious differences between the calls

You could provide intensities from different SNPs too though I dont' think there's a point.


# idat2cluster.py

IDAT files contain the red/green intensities (and other data) from a genotyping run. Each IDAT file contains the calls of one individual. This code does two things: first it extracts out the red/green intensities into human readable text  files (that can be processed by _plotcluster.py` above); and second the output is reorganised so that instead of being organsied by person, they are organised by SNP. That is the input has one file per individual which each contains all the SNP data, while the output has one file per SNP which each contains all the individual data

### Requirements

python3 and numpy

### Inputs

The following inputs are required
* Optionally a sample sheet (names of individuals with the mapping to the IDs in the IDAT files). If you don't have it, put `NONE`
* A  directory that contains  IDAT files. If there are nested directories all nested directories will be included too
* A chip manifest file in CSV format 
* An output directory, specified with the `-o` option. This is compulsory.
* number of threads to use `-n` -- optional. This has not been benchmarked but modest performance improvement can be gained by choosing up to 8.
* There are other options that will be used at some point


### Outputs

The code produces millions of output files, one for each SNP on the chip. Putting this all in one directory would be unmanageable. The output directory is organised hierarchically as follows:
* One directory for each chromosome
   * within that directory, SNPs are divided into bins by chromosome position. The bins are 100k in size. So `output/12/1203` consists lf all the SNPs that are found on chromosome 12 between positions $1203\times 100,000$ and $1204\times 100,000-1$ that is between positions 120,300,000 and 120,399,999.
        * within each of these directories you will find CSV files named by SNP (names given in the chip manifest)

You will need a manifest file to navigate this.

