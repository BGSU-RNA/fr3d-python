FROM python:2.7

WORKDIR /usr/src/app

RUN \
    wget https://mmcif.wwpdb.org/docs/sw-examples/python/src/pdbx.tar.gz && \
    tar -xzf pdbx.tar.gz && \
    rm pdbx.tar.gz && \
    mkdir -p source/python/modules && \
    mv pdbx source/python/modules

ENV PYTHONPATH "${PYTHONPATH}:/usr/src/app/source/python/modules"

WORKDIR /rna

ADD . /rna

RUN pip install -r requirements-python-2-7.txt

CMD ["/bin/bash"]
