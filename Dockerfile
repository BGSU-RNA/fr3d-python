FROM python:3.9

WORKDIR /rna

ADD . /rna

RUN python setup.py install

ENTRYPOINT ["/bin/bash"]
