FROM continuumio/miniconda3

WORKDIR /usr/src/caddie/

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt-get update
RUN apt-get -y install build-essential
RUN apt-get -y install gcc

RUN conda install -c conda-forge -y graph-tool

RUN apt-get update
RUN apt-get install -y supervisor nginx
RUN apt-get install wget

COPY ./requirements.txt /usr/src/caddie/requirements.txt

# RUN pip3 install setuptools==57
RUN pip3 install -r /usr/src/caddie/requirements.txt
RUN pip3 install gunicorn

COPY ./supervisord.conf /etc/supervisor/conf.d/supervisord.conf
COPY ./docker-entrypoint.sh /entrypoint.sh
COPY ./import-data.sh /import-data.sh

COPY . /usr/src/caddie/

EXPOSE 8000

ENTRYPOINT ["sh", "/entrypoint.sh"]
