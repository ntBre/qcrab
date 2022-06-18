TESTFLAGS = --nocapture --test-threads=1

test:
	cargo test -- ${TESTFLAGS}
