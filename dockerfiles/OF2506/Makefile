ME=$(shell whoami)
clean:
	rm -f hello

build:
	docker build -t ithacafv/openfoam2506-muq2-pytorch -f ./Dockerfile .

run: 
	docker run -it ithacafv/openfoam2506-muq2-pytorch bash