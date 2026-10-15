FROM python:3.11

WORKDIR /rna

ADD . /rna

RUN python -m pip install .
RUN pip install biopython

CMD ["/bin/bash"]
