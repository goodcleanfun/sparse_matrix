
test:
	clib install --dev --concurrency 1
	@$(CC) test.c -std=c99 -I src -I deps -o $@
	@./$@

.PHONY: test
