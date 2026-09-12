.PHONY: test-python test-cpp test-web test build-cpp build-web

test-python:
	cd python && (test -d .venv || python3 -m venv .venv) && \
	  .venv/bin/pip install -q -r requirements.txt && \
	  .venv/bin/python -m pytest -q

build-cpp:
	cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release
	cmake --build cpp/build -j

test-cpp: build-cpp
	./cpp/build/wavefield_cpp --nx 32 --ny 32 --steps 20

test-web:
	cd web && npm install && npm test && npm run build

test: test-python test-cpp test-web
