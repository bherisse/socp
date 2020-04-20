COVID-19
===========

covid-19 model is in src/models/covid19
main program is in tests/testCovid19.cpp

Using a SEIR model for the covid-19 epidemic, the problem considered here is the following:
- Control lockdown with minimum effort (in a quadratic way).
- One year after end of lockdown (2020, May 11 in France), the total removed (recovered + isolated + fatalities) population achieves 80% (which corresponds approximatively to natural herd immunity).

Even for this "worst case" scenario (80% herd immunity will induce many fatalities), results confirm that end of lockdown sould be very progressive to restrain the number of infections actively circulating along the next year. 

Note that initial covid-19 parameters should be refined.

Future work: control the reproduction number to minimize lockdown effort while constraining the number of infected people.