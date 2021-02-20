perform_tests:
	R -e "library('devtools')"
	R -e "devtools::test()"

perform_cran_check:
	R -e "library('devtools')"
	R -e "devtools::check()"
