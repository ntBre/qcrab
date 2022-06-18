TESTFLAGS = --nocapture

test:
	cargo test -- ${TESTFLAGS}
