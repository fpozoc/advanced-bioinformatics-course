# Software requirements

Remember that all the packages required to run the analysis are already installed in our [Docker image](https://hub.docker.com/r/osvaldogc/ufv). 

Install the following packages if you want to run your analysis locally.

```bash
#! /bin/bash
SOFTWARE=$PWD/software
mkdir -p $SOFTWARE 
echo "$SOFTWARE created!"
```

## Install FastQC 0.11.5

```bash
echo "Downloading FastQC..."
wget -q -P $SOFTWARE/ http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip -q $SOFTWARE/fastqc_v0.11.5.zip -d $SOFTWARE && rm $SOFTWARE/fastqc*.zip
echo "FastQC downloaded here: $SOFTWARE"
FASTQC=$SOFTWARE/FastQC
export PATH=$FASTQC:$PATH
chmod 755 $FASTQC
echo "fastqc executable has been located here: $FASTQC"
```

## Install Trimmomatic 0.36

```bash
echo "Downloading Trimmomatic..."
wget -q -P $SOFTWARE/ http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
echo "Trimmomatic downloaded here: $SOFTWARE"
unzip -q $SOFTWARE/Trimmomatic-0.36.zip -d $SOFTWARE && rm $SOFTWARE/Trimmomatic-0.36.zip
echo "Trimmomatic is ready to be executed with: java -jar $SOFTWARE/Trimmomatic-0.36/trimmomatic-0.36.jar --help"
```

## Install Bowtie 2

```bash
echo "Downloading Bowtie 2..."
wget -q -P $SOFTWARE/ https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip
echo "Bowtie 2 downloaded here: $SOFTWARE"
unzip -q $SOFTWARE/bowtie2-2.3.3.1-linux-x86_64.zip -d $SOFTWARE && rm $SOFTWARE/bowtie2-2.3.3.1-linux-x86_64.zip
BOWTIE2=$SOFTWARE/bowtie2-2.3.3.1-linux-x86_64
export PATH=$BOWTIE2:$PATH
echo "Bowtie 2 executable has been located here: $FASTQC"
```

## Install Samtools 1.9

```bash
wget -q -P $SOFTWARE/ https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar xjf $SOFTWARE/samtools-1.9.tar.bz2 -C $SOFTWARE && rm $SOFTWARE/samtools-1.9.tar.bz2
cd $SOFTWARE/samtools-1.9/ && ./configure --prefix=$SOFTWARE/samtools-1.9
cd $SOFTWARE/samtools-1.9/ && make
cd $SOFTWARE/samtools-1.9/ && make install
SAMTOOLS=$SOFTWARE/samtools-1.9/bin
export PATH=$SAMTOOLS:$PATH
```

## Install Stringtie 1.3

```bash
echo "Downloading Stringtie..."
wget -q -P $SOFTWARE/ http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.0.Linux_x86_64.tar.gz
echo "Stringtie downloaded here: $SOFTWARE"
tar xvfz $SOFTWARE/stringtie-1.3.0.Linux_x86_64.tar.gz -C $SOFTWARE && rm $SOFTWARE/stringtie-1.3.0.Linux_x86_64.tar.gz
STRINGTIE=$SOFTWARE/stringtie-1.3.0.Linux_x86_64
export PATH=$STRINGTIE:$PATH
echo "stringtie executable has been located here: $STRINGTIE"
```
