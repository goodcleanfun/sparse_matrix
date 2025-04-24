
test:
	clib install --dev --concurrency 1
	@$(CC) test.c $(CFLAGS) -I src -I deps -o $@
	@./$@

.PHONY: test
