VERSION = 1.0

.PHONY: build
build:
	go mod download
	go build -o bin/depleteIndiv     ./cmd/depleteIndiv
	go build -o bin/depleteMarkers   ./cmd/depleteMarkers
	go build -o bin/importHaploCall  ./cmd/importHaploCall
	go build -o bin/markerOverlap    ./cmd/markerOverlap
	go build -o bin/pseudoHaploidize ./cmd/pseudoHaploidize
	go build -o bin/filterRelates    ./cmd/filterRelates
	go build -o bin/simContam        ./cmd/simContam
