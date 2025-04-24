
test:
	clib install --dev --concurrency 1
	@$(CC) test.c $(CFLAGS) -I src -I deps $(LDFLAGS) -o $@
	@./$@

.PHONY: test
