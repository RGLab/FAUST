perform_unit_tests:
	R -e "library('devtools')"
	R -e "devtools::check()"
