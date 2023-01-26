# syntax=docker/dockerfile:1
# FROM python:3.8-slim-bullseye as base
FROM python:3.8-alpine3.16 as base

WORKDIR /

ENV PYTHONUNBUFFERED: 1

RUN uname -r

# RUN apt-get update 
# # RUN apt-get dist-upgrade
# # RUN uname -r

# RUN apt-get install -y build-essential
# RUN apt-get install -y libz-dev
# RUN apt-get install -y bedtools

RUN apk up-date
RUN apk upgrade
RUN apk add build-essential
RUN apk add libz-dev
RUN apk add bedtools


COPY requirements.txt requirements.txt
COPY sge-primer-scoring/requirements.txt scoring_requirements.txt

RUN pip3 install -U pip wheel setuptools 
RUN pip3 install -r requirements.txt
RUN pip3 install -r scoring_requirements.txt

COPY . .

RUN df -h /tmp

FROM base as dev
ENTRYPOINT [ "python3", "./src/cli.py" ]

FROM base as test
CMD ["python3", "-m", "unittest"]


# TODO: Flesh out loading data into Docker Container ticket
