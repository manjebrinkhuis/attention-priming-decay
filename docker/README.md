# neuro-docker

Dockerized environment to use `nipype`, `fsl` and
`freesurfer` (as a mounted volume) and jupyterlab.

On Linux, `./run_terminal.sh` runs a container, with a mounted
data directory (change the `/path/to/your/data` according to your system),
while sharing your X11 socket with docker, so you can also run graphical
applications `freeview` and `fslview`.

## Usage

A workflow could be to run a container, run `jupyter lab` in the container
and browse to `localhost:8888`, to work in jupyterlab.

Note that you can also run `fslview` and `freeview` from a terminal
in jupyterlab.
