FROM ubuntu:16.04
MAINTAINER Ewa Piskadlo <ewa.piskadlo@protonmail.com>


RUN apt-get update 
RUN apt-get install -y \
	python2.7 \
	python-pip 

RUN pip install --upgrade pip

RUN pip install biopython==1.70

COPY blockParse.py . 
	
