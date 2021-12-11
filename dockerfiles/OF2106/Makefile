ME=$(shell whoami)
clean:
	rm -f hello

build:
	docker build --build-arg myuser=$(shell whoami) --build-arg UID=$(id -u) --build-arg GID=$(id -g) -t ithacafv/openfoam2106 -f ./Dockerfile .

run: build
	xhost +"local:docker" && \
	docker run -ti --rm \
        -v "${HOME}:/home/${ME}" \
		-e DISPLAY=${DISPLAY} \
        --security-opt seccomp=unconfined \
		ubuntu-apc


