# How to Update the Official Docker Image

To update the official Docker image on `registry.dune-project.org`,
please follow these steps:

1. Actually build the docker image:
  ```
  docker build . --no-cache -t registry.dune-project.org/dune-fem/dune-fempy
  ```
  Please do not omit `--no-cache` to ensure every step of the build is actually
  executed.

  *Note*: This will set up the Debian base system of the Docker image, which takes
  quite some time.

1. Verify your build by testing the Docker image locally:
   ```
   docker run --rm -v dune:/dune -p 127.0.0.1:8888:8888 \
      registry.dune-project.org/dune-fem/dune-fempy
   ```
   Use your favorite web browser, visit http://127.0.0.1:8888, enter the password
   `dune`, and run the demo notebooks (depending on the jupyter version you
   might be given  tokenized http address to use instead).
   *Note*: If you already have a volume `dune`, you will not receive the updated demo notebooks. It might be better to start with a clean volume.

1. Log into the Docker registry:
  ```
  docker login registry.dune-project.org
  ```

1. Push the Docker image:
  ```
  docker push registry.dune-project.org/dune-fem/dune-fempy
  ```
  This step actually updates the Docker image on the server.

  **Warning:** Updating the official Docker image affects all users of the
  image once they pull it.
  Make sure you don't push a buggy combination of the DUNE modules.

1. Log out of the Docker registry:
  ```
  docker logout registry.dune-project.org
  ```

That's it.
