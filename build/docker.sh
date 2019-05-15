docker build\
  --build-arg userId=$(id -u) --build-arg groupId=$(id -g)\
  -t dune .

docker run -it -v $HOST:/host -v dune:/dune\
  -v /tmp/.X11-unix:/tmp/.X11-unix:ro -e DISPLAY=:0\
  dune bash
