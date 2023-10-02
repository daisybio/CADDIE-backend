FROM continuumio/miniconda3:4.10.3

WORKDIR /usr/src/caddie/

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# install gcc
RUN apt-get update
RUN apt-get -y install build-essential
RUN apt-get -y install python-dev python3-dev
RUN apt-get -y install gcc

# RUN apt-get update
# RUN apt-get install -y supervisor nginx
# RUN apt-get install -y libgtk-3-dev
# RUN apt-get install -y wget


RUN conda install -y conda=4.10.3
RUN conda install -c conda-forge -y graph-tool

# setuptools 58 does not support use_2to3
# RUN pip install setuptools==57.5.0
# RUN pip install -r /usr/src/caddie/requirements.txt
# RUN pip install gunicorn




# OLD
RUN apt-get update
RUN apt-get install -y supervisor nginx
RUN apt-get install -y libgtk-3-dev
RUN apt-get install wget

COPY ./requirements.txt /usr/src/caddie/requirements.txt

# RUN conda install -y conda=4.3.16
# RUN conda install -c conda-forge -y graph-tool=2.32

RUN pip install -r /usr/src/caddie/requirements.txt
RUN pip install gunicorn






# RUN pyensembl install --release 106 --species human

COPY ./supervisord.conf /etc/supervisor/conf.d/supervisord.conf
COPY ./docker-entrypoint.sh /entrypoint.sh
COPY ./import-data.sh /import-data.sh

COPY . /usr/src/caddie/

EXPOSE 8000

ENTRYPOINT ["sh", "/entrypoint.sh"]
