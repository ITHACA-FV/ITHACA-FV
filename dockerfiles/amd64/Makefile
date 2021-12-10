ME=$(shell whoami)

build:
	docker build --no-cache  -t ithacafv/ithacafv:manifest-amd64 -f ./Dockerfile .

run:
	docker run -ti --rm \
		-v "${HOME}:/home/ithacafv/${ME}" \
        --security-opt seccomp=unconfined \
		ithacafv/ithacafv:manifest-amd64


