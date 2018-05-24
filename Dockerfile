###############################################
# Dockerfile to build readMapper container image
# Based on Ubuntu 16.04
# Build with:
#   sudo docker build -t readmapper .
###############################################

    # Use ubuntu 16.04 base image
    FROM ubuntu:16.04

    # File author/maintainer info
    MAINTAINER Aurélien BIRER <abirer@chu-clermontferrand.fr>

    # set non-interactive mode
    ENV DEBIAN_FRONTEND noninteractive

    #From ariba docker
    RUN apt-get update
    RUN apt-get install --no-install-recommends -y \
        build-essential \
        cd-hit \
        curl \ 
        git \ 
        libbz2-dev \
        liblzma-dev \
        mummer \
        python \ 
        python-pip \
        python-dev \
        python-setuptools \
        python3-dev \
        python3-setuptools \
        python3-pip \
        python3-tk \
        python3-matplotlib \
        unzip \
        wget \
        zlib1g-dev

    #For Resolve encoding problems with docker
    RUN apt-get update --fix-missing 
    RUN apt-get install --no-install-recommends -y locales

    # Set UTF8 for system
    RUN locale-gen en_US.UTF-8
    ENV LANG en_US.UTF-8
    ENV LC_ALL en_US.UTF-8

    #resolve pip3 problem
    RUN    curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
           python3 get-pip.py --force-reinstall


    #Install readmapper dependencies
    RUN     pip install --upgrade pip && \
	    pip install pandas && \
            pip install openpyxl && \
            pip install biopython && \
            pip install python-docx && \
            pip install xlrd && \
	    pip install Pillow


    RUN    wget -q http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip &&\
           unzip bowtie2-2.2.9-linux-x86_64.zip  &&\
           rm bowtie2-2.2.9-linux-x86_64.zip 

    # Need MPLBACKEND="agg" to make matplotlib work without X11, otherwise get the error 
    # _tkinter.TclError: no display name and no $DISPLAY environment variable 
    ENV ARIBA_BOWTIE2=$PWD/bowtie2-2.2.9/bowtie2 ARIBA_CDHIT=cdhit-est MPLBACKEND="agg" 
    

    RUN git clone https://github.com/sanger-pathogens/ariba.git  &&\
        cd ariba  &&\
        git checkout v2.10.1  &&\
        python3 setup.py test  &&\
        python3 setup.py install

    CMD ariba

    RUN    cd /usr/local  &&\
	   git clone https://github.com/CNRResistanceAntibiotic/readMapper.git && \
	   python3 readMapper/setup.py install
    
    
    RUN    apt autoremove --purge --yes && \
           apt clean
