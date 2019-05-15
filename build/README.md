Docker
======
To build the container run
```
docker build --build-arg userId=$(id -u) --build-arg groupId=$(id -g) -t dune .
```
and to get into the container with the current directory are working
directory execute
```
docker run -it -v $PWD:/host -v dune:/dune\
  -v /tmp/.X11-unix:/tmp/.X11-unix:ro -e DISPLAY=:0\
  dune bash
```

Vagrant
=======
Run `vagrant up` to get started and then `vagrant ssh` in your working
directory to go into the container.
To update the container run `vagrant provision` which will reexecute the
`bootstrap.sh` file.

Suggestions
===========
Execute the `build` dune script in the bootstrap script putting the dune
folders into the `home` directory - set the `DUNE_CONTROL_PATH` to
`$HOME:.` in `bashrc`. Also define some commands to update the `dune-py`
(i.e. call `dune-setup.py`) and also to generate a new virtual env with
it's own `dune-py` folder.
