TESTFLAGS = --nocapture --test-threads=1
ARGS =

test:
	cargo test -- ${TESTFLAGS} $(ARGS)

clippy:
	cargo clippy --tests
