############################################################
# Singularity to build readMapper pipeline
# Based on Ubuntu 16.04
#   sudo singularity build readMapper.simg Singularity 
############################################################

# Set the base format to docker
Bootstrap: docker

# Set the base image to Ubuntu
From: ubuntu:16.04

%help
    Help me. I'm in the container.

%labels
    Maintainer Aurélien BIRER <abirer@chu-clermontferrand.fr>
    Version v1.0

%post

    #From ariba docker
    apt-get update
    apt-get install --no-install-recommends -y build-essential cd-hit curl git libbz2-dev liblzma-dev mummer python python-pip python-dev python-setuptools python3-dev python3-setuptools python3-pip python3-tk python3-matplotlib unzip wget zlib1g-dev libcurl4-openssl-dev

    #For Resolve encoding problems with docker
    apt-get update --fix-missing 
    apt-get install --no-install-recommends -y locales

    # Set UTF8 for system
    locale-gen en_US.UTF-8

    #resolve pip3 problem
    curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
        python3 get-pip.py --force-reinstall

    #Install readmapper dependencies
    pip install --upgrade pip && \
        pip install pandas && \
        pip install openpyxl && \
        pip install biopython && \
        pip install python-docx && \
        pip install xlrd && \
        pip install Pillow

    wget -q http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip &&\
        unzip bowtie2-2.2.9-linux-x86_64.zip  &&\
        rm bowtie2-2.2.9-linux-x86_64.zip 

    #ariba env

    # Need MPLBACKEND="agg" to make matplotlib work without X11, otherwise get the error 
    # _tkinter.TclError: no display name and no $DISPLAY environment variable

    echo 'export ARIBA_BOWTIE2=/bowtie2-2.2.9/bowtie2' >>$SINGULARITY_ENVIRONMENT
    export ARIBA_BOWTIE2=/bowtie2-2.2.9/bowtie2
    echo 'export MPLBACKEND="agg"' >>$SINGULARITY_ENVIRONMENT
    export MPLBACKEND="agg"
    echo 'export ARIBA_CDHIT=cdhit-est' >>$SINGULARITY_ENVIRONMENT
    export ARIBA_CDHIT=cdhit-est

    git clone https://github.com/sanger-pathogens/ariba.git  &&\
        cd ariba  &&\
        git checkout v2.10.1  &&\
        python3 setup.py test  &&\
        python3 setup.py install

    cd /usr/local  &&\
        git clone https://github.com/CNRResistanceAntibiotic/readMapper.git && \
        cd readMapper && \
        python3 setup.py install

    apt autoremove --purge --yes && \
        apt clean


%environment

    # Set UTF8 for system
    LANG=en_US.UTF-8
    LC_ALL=en_US.UTF-8


