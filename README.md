# spots2reciprocal
Script to transform a list of 2D Bragg rflection centroids to the reciprocal space.

## Overview 

This script is meant to processes information from 2D Bragg's reflections centroids from from [XDS](https://xds.mr.mpg.de/html_doc/xds_files.html) (i.e. [XDS_ASCII.HKL](https://xds.mr.mpg.de/html_doc/xds_files.html#XDS_ASCII.HKL) and [SPOT.XDS](https://xds.mr.mpg.de/html_doc/xds_files.html#SPOT.XDS)) and transform them to 3D using [spot2pdb](https://strucbio.biologie.uni-konstanz.de/pub/linux_bin/spot2pdb).

## Installation

Linux is fully supported. 

Download [spot2pdb](https://strucbio.biologie.uni-konstanz.de/pub/linux_bin/spot2pdb) and put it in some path where it's globally executable.

```bash
wget https://strucbio.biologie.uni-konstanz.de/pub/linux_bin/spot2pdb
mv spot2pdb /usr/local/bin
sudo chmod a+x /usr/local/bin/spot2pdb
```

Create a conda environment using the provided [environment.yml](./environment.yml) file:

```bash
conda env create -n spots2reciprocal --file environment.yml
```

## Usage

Run it with...