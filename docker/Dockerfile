FROM debian:stretch
MAINTAINER "Martin Nolte"

RUN apt-get update && apt-get dist-upgrade --yes --no-install-recommends
RUN apt-get install --yes --no-install-recommends \
  ca-certificates wget tar \
  libc6-dev gcc g++ binutils make cmake pkg-config patch \
  python3-dev python3-setuptools python3-pip python3-wheel \
  python3-pandocfilters \
  python3-numpy python3-scipy python3-matplotlib \
  python3-ufl

RUN pip3 install jupyter

ADD install-dune.sh documentation.patch /root/
RUN /root/install-dune.sh

RUN adduser --disabled-password --home /home/dune --gecos "" --uid 50000 dune
USER dune
WORKDIR /home/dune

RUN cp -r /usr/local/share/doc/dune-corepy/jupyter /home/dune/dune-corepy
RUN cp -r /usr/local/share/doc/dune-fempy/jupyter /home/dune/dune-fempy

ENV LD_LIBRARY_PATH "/usr/lib:/usr/local/lib"
ENV DUNE_CMAKE_FLAGS="CMAKE_BUILD_TYPE=Release BUILD_SHARED_LIBS=TRUE DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS=TRUE"

EXPOSE 8888
CMD ["/usr/local/bin/jupyter-notebook", "--no-browser", "--ip=*"]
