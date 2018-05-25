# Create conda volume
# docker volume create neuro

#
docker run --rm \
    -it \
    -p 8888:8888 \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -e DISPLAY=unix$DISPLAY \
    -v /mnt/ext_data:/data/ext \
    -v /mnt/data:/data/host \
    -v /mnt/data/Apps/Linux/freesurfer:/opt/freesurfer \
    --name neuro-docker \
    manjebrinkhuis/neuro-docker
