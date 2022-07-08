# Course in Advanced Bioinformatics

## Table of Contents

* [Outline](#outline)
* [Software requirements](#software-requirements)
    * [Running with our Docker Image](#running-with-our-docker-image)
    * [Running locally](#running-locally)
* [Some EXTRA recommendations](#some-extra-recommendations)
* [Author](#author)
* [Contributing](#contributing)

# Outline

|#  |Date    |Topic                                                                                     |Lecturer     |Theory                                                                                                                                                                                                                                                                                                |Hands-on                                                                                                                                                                                                     |
|---|--------|------------------------------------------------------------------------------------------|-------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|1  |01/09/21|Introduction to genome sequencing and next generation sequencing methods                  |Fernando Pozo|[pdf](https://gitlab.com/fpozoc/advanced_bioinformatics/-/blob/master/ADB-01-int-bioinf/theory/main.pdf)                                                                                                                                                                                              |-                                                                                                                                                                                                            |
|2  |02/09/21|Introduction to the Linux command line, NGS data formats, read mapping and alignments     |Fernando Pozo|[pdf](https://gitlab.com/fpozoc/advanced_bioinformatics/-/blob/master/ADB-02-ngs-linux/theory/main.pdf)/[web](https://gitlab.com/fpozoc/advanced_bioinformatics/-/tree/master/ADB-02-ngs-linux/theory)                                                                                                |[pdf](https://gitlab.com/fpozoc/advanced_bioinformatics/-/tree/master/ADB-02-ngs-linux/handson/hands-on.pdf)/[web](https://gitlab.com/fpozoc/advanced_bioinformatics/-/tree/master/ADB-02-ngs-linux/handson)|
|3  |10/09/21|Introduction to Genome Browsing with the UCSC Genome Browser                              |Osvaldo Graña|[pdf](https://gitlab.com/fpozoc/advanced_bioinformatics/-/blob/master/ADB-03-gen-brow/theory/main.pdf)                                                                                                                                                                                                |[pdf1](https://gitlab.com/fpozoc/advanced_bioinformatics/-/blob/master/ADB-03-gen-brow/hands-on/hands-on-1.pdf)/[pdf2](https://gitlab.com/fpozoc/advanced_bioinformatics/-/blob/master/ADB-03-gen-brow/hands-on/hands-on-2.pdf)|
|4  |09/09/21|NGS applications: Sequence assembly, de‐novo sequencing and EST sequencing                |Osvaldo Graña|[pdf](https://gitlab.com/fpozoc/advanced_bioinformatics/-/blob/master/ADB-04-seq-asse/theory/main.pdf)                                                                                                                                                                                                |-                                                                                                                                                                                                            |
|5  |08/09/21|NGS applications: ChIP‐seq                                                                |Osvaldo Graña|[pdf](https://gitlab.com/fpozoc/advanced_bioinformatics/-/blob/master/ADB-05-chip-seq/theory/main.pdf)                                                                                                                                                                                                |[pdf](https://gitlab.com/fpozoc/advanced_bioinformatics/-/blob/master/ADB-05-chip-seq/hands-on/hands-on.pdf)                                                                                                 |
|6  |03/09/21|NGS applications: Variant detection, SNPs, CNVs, structural variants, epigenetic variation|Elena Piñeiro|[pdf](https://gitlab.com/fpozoc/advanced_bioinformatics/-/blob/master/ADB-06-variants/theory/main.pdf)                                                                                                                                                                                                |-                                                                                           |
|7  |06/09/21|NGS applications: RNA‐seq                                                                 |Fernando Pozo|[web](https://gitlab.com/fpozoc/advanced_bioinformatics/-/tree/master/ADB-07-rna-seq/theory)                                                                                                                                                                                                          |[web](https://gitlab.com/fpozoc/advanced_bioinformatics/-/tree/master/ADB-07-rna-seq/hands-on)                                                                                                               |


# Software requirements

## Running with our Docker Image

* [osvaldogc/ufv:2.0](https://hub.docker.com/r/osvaldogc/ufv/): Docker image with `Hands On` packages already installed.
Images can be pulled via `docker pull osvaldogc/ufv:2.0`. Detailed instructions [here](https://gitlab.com/fpozoc/advanced_bioinformatics/-/tree/master/ADB-02-ngs-linux/handson#21-first-steps-with-docker).

## Running locally

Remember that all the packages required to run the analysis are already installed in our [Docker image](https://hub.docker.com/r/osvaldogc/ufv). 

Follow this [guide](https://gitlab.com/fpozoc/advanced_bioinformatics/-/tree/master/ADB-02-ngs-linux/handson/local_installation) if you want to run your analysis locally.

# Some EXTRA recommendations

If you want to complement, add or step up your knowledge in bioinformatics... 

Some deeper and longer... but highly recommended bioinformatics courses (still available online):

- [Foundations of Computational and Systems Biology](https://ocw.mit.edu/courses/biology/7-91j-foundations-of-computational-and-systems-biology-spring-2014/). MIT Opencourseware. Spring 2014.
- [MIT CompBio](https://www.youtube.com/watch?v=sX4cMu9Azgs&list=PLypiXJdtIca6GBQwDTo4bIEDV8F4RcAgt). MIT. 2018.

Recommended books:

- [Bioinformatics algorithms](https://www.bioinformaticsalgorithms.org/read-the-book) by Pavel A. Pevzner.
- [Modern Statistics for Modern Biology](https://www.huber.embl.de/msmb/) by Susan Holmes and Wolfgang Huber.

# Author

Fernando Pozo – [@fpozoca](https://twitter.com/fpozoca) - [Google Scholar](https://scholar.google.com/citations?user=3YLw4PQAAAAJ&hl=en) – fpozoc@cnio.es - [gitlab.com/fpozoc](https://gitlab.com/fpozoc)

Thanks to Gonzalo Gomez and Coral Fustero for providing me highly interesting notes and materials for the Introduction to NGS and RNA-seq class, and for the rest of the [Bioinformatics Unit](https://bioinformatics.cnio.es/staff/fpozoc/) ([@BU_CNIO](https://twitter.com/BU_CNIO)) staff which is based on 
[Spanish National Cancer Research Centre (CNIO)](https://www.cnio.es/).


# Contributing

1. Fork it (<https://github.com/fpozoc/advanced_bioinformatics>).
2. Create your feature branch (`git checkout -b feature/fooBar`).
3. Commit your changes (`git commit -am 'Add some fooBar'`).
4. Push to the branch (`git push origin feature/fooBar`).
5. Create a new Pull Request.
