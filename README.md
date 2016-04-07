# made: Microarray Analysis of Differential Expression

This package has completed development and is being submitted to Bioconductor.

In the interim, download the latest version from the 
[releases](https://github.com/fboulnois/made/releases) tab and extract the zip
file to a folder.

**Note:** Before starting the installation of `made`, you may need to install
`libxml2-dev` or equivalent from the terminal if using a Linux distribution:

- Type `sudo apt-get install libxml2-dev` in Debian or Ubuntu
- Type `yum install libxml2-devel` in Fedora or Red Hat Enterprise Linux
- Type `sudo zypper install libxml2-devel` in OpenSUSE

Then, in R:

1. Use `setwd()` to set the working directory to that folder
2. Run `source("madeinstall.R")` to automatically install the required
dependencies from CRAN and Bioconductor.

If instead you are using RStudio then you can open `madeinstall.R` and source it
 directly.

**Note:** ReactomePA and GOstats are large dependencies and may take a while to 
download. Annotation packages must also be downloaded to support a particular 
microarray platform.
