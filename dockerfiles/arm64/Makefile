ME=$(shell whoami)

build:
	docker build --build-arg myuser=$(shell whoami) --build-arg UID=$(id -u) --build-arg GID=$(id -g) -t ithacafv/ithacafv:manifest-arm64 -f ./Dockerfile .

run:
	docker run -ti --rm \
		-v "${HOME}:/home/ithacafv/${ME}" \
        --security-opt seccomp=unconfined \
		ithacafv/ithacafv:manifest-arm64


