# Instruction to build and push docker image

To build image, execute command

    $ docker build .

in the current directory. `docker` will read file `./Dockerfile` and build an
image from it.
Sometimes, there can be cache problems, when the image is build with outdated
cache data. In such cases, flag `--no-cache` can be used:

    $ docker build --no-cache .


## Tagging an image

After successful build, the image must be tagged (or versioned)

	docker tag <image-id> swagkaust/swagimage:YYYY-MM.BN

where <image-id> is ID of the built image, which is provided in the output of
the `docker build` command (in the last line), `YYYY` is current year,
`MM` the current month with leading zero, `BN` the build number.
For example, the very first build in August 2018 was 2018-08.01.

It is a common practice to also put a tag `latest` on the most recent image:

    docker tag <image-id> swagkaust/swagimage:latest 


## Upload Docker image to the cloud

To push image to the Docker cloud (registry), execute command

    docker push swagkaust/swagimage:latest
    docker push swagkaust/swagimage:YYYY-MM.BN

The repository is located at:
https://cloud.docker.com/swarm/swagkaust/repository/list

