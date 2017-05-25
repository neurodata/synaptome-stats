# Commencing 
FROM ubuntu:16.04

## Updates and software that I like
# Install vnc, xvfb in order to create a 'fake' display and firefox
RUN apt-get update  && apt-get install -y \
  man \
  ssh \
  curl \
  python

RUN apt-get install -y python-pip python-dev build-essential
RUN pip install --upgrade pip

RUN pip install numpy
RUN pip install ndio==1.1.5
RUN pip install scipy
RUN pip install h5py

#### Add user
RUN useradd -ms /bin/bash py
RUN chown py:py /home/py & addgroup py staff

RUN mkdir /home/py/scripts
RUN mkdir /home/py/data
RUN mkdir /home/py/output

COPY ./getF0.py /bin/getF0.py
COPY ./functions.py /bin/functions.py
COPY ./runFeatureF0.sh /bin/runFeatureF0.sh
COPY ./runFeatureF0.sh /runFeatureF0.sh

RUN chown py:py /bin/getF0.py /bin/functions.py
RUN chown py:py -R /home/py/data
RUN chown py:py -R /home/py/scripts
RUN chown py:py -R /home/py/output

RUN chmod +x /bin/getF0.py
RUN chmod +x /bin/functions.py
RUN chmod +x /bin/runFeatureF0.sh

USER py
WORKDIR /home/py/

ENTRYPOINT ["runFeatureF0.sh"]

