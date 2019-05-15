FROM ubuntu:18.04
MAINTAINER "Andreas Dedner"

ARG userId
ARG groupId
RUN addgroup --gid $groupId dune
RUN adduser --disabled-password --home /dune --gecos "" --uid $userId --gid $groupId dune

COPY bootstrap.sh /
RUN bash /bootstrap.sh

COPY buildDUNEFemDG.sh /dune/
RUN chown dune:dune /dune/buildDUNEFemDG.sh
RUN chmod 770 /dune/buildDUNEFemDG.sh

USER dune
WORKDIR /dune
RUN ./buildDUNEFemDG.sh

ENV DISPLAY :0
EXPOSE 8888
VOLUME ["/dune"]
WORKDIR /host
