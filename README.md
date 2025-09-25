## containers
For each tool in the repository, the corresponding Docker container was build from scratch starting from a Dockerfile in order to pursuit version control (both sfotware and its dependencies) and data reproducibility.

### to build a container from a Dockerfile
The less the number of layers used to build the container, the better. 

```bash
# login into your account
docker login -u <username> -p <password> docker.io

# build the container - you must be in the same directory as the Dockerfile
docker build -t username/tool:version .
```

### basic commands to handle a Docker container

```bash
# check inside the container
docker run --rm -v path-to-wd/:/tmp --entrypoint "/bin/bash" -it username/tool:version

# push the container into the hub
docker oush username/tool:version

# pull a docker container - the container must be public
docker pull username/tool:version

# check all images
docker images

# remove unused images
docker rmi <tool_tag>
```
### reference
[Docker Hub.](https://hub.docker.com/)
