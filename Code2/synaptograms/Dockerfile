# Commencing 
FROM neurodata/meda:v0.1.2.9009

## Updates and software that I like
# Install vnc, xvfb in order to create a 'fake' display and firefox
USER root
RUN apt-get update  && apt-get install -y \
  man \
  ssh \
  curl \
  python

RUN apt-get install -y python-pip python-dev build-essential
RUN apt-get install -y python3-pip python3-dev
RUN pip install --upgrade pip

RUN pip3 install numpy
RUN pip3 install ndio==1.1.5
RUN pip3 install scipy
RUN pip3 install h5py

COPY ./runSynaptograms.sh /runSynaptograms.sh

WORKDIR /home/meda
RUN mkdir scripts
RUN mkdir data

COPY ./functions.py /bin/functions.py
COPY ./tryCubes.py /bin/tryCubes.py
COPY ./k15Synaptograms.R /bin/k15Synaptograms.R
COPY ./runSynaptograms.sh /bin/runSynaptograms.sh

RUN chown meda:meda -R /home/meda/data
RUN chown meda:meda -R /home/meda/scripts

RUN chmod +x /bin/runSynaptograms.sh
RUN chmod +x /bin/runSynaptograms.sh
RUN chmod +x /bin/tryCubes.py
RUN chmod +x /bin/k15Synaptograms.R 
RUN chmod +x /bin/functions.py

USER meda

ENTRYPOINT ["runSynaptograms.sh"]

