ME=$(shell whoami)
clean:
	rm -f hello

build:
	docker build -t ithacafv/openfoam2412-muq2-pytorch -f ./Dockerfile .

run: build
	xhost +"local:docker" && \
	docker run -ti --rm \
        -v "${HOME}:/home/${ME}" \
		-e DISPLAY=${DISPLAY} \
        --security-opt seccomp=unconfined \
		ubuntu-apc


