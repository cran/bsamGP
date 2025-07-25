
# bsamGP README

## Authors

	* Seongil Jo <statjs@jbnu.ac.kr>
	* Taeryon Choi <trchoi@korea.ac.kr>
	* Beomjo Park <beomjop@alumni.cmu.edu>
	* Peter J. Lenk <plenk@umich.edu>

## Acknowledgments

	* Research of Seongil Jo was supported by Basic Science Research Program
	through the National Research Foundation of Korea (NRF) funded by
	the Ministry of Education (NRF-2017R1D1A3B03035235)

	* Research of Taeryon Choi was supported by Basic Science Research Program
	through the National Research Foundation of Korea (NRF) funded by
	the Ministry of Education (NRF-2016R1D1A1B03932178)


## Changelog
	* v1.0.1 -- 1.0.2
		* Add missing dependencies to src/Makevars. (Thanks to Prof. Brian Ripley)

	* v1.1.0
		* Change input UI to R formula.
		* Add S3methods : predict, summary
		* Clean plot method.

	* v1.1.1 -- 1.1.3
		* Separate predictbsam & predictgbsam : BUGFIX
		* Add verbose option.
		* Minor Bug fix
		
	* v.1.2.0 -- 1.2.7
		* Supports Multiple extreme shapes for bsar() and bsaq()
		* Supports scalable regression function bsarBig()
		* Memory bug fix
    * Minor bug fix
    * gcc/gfortran 12 compatibility check and gcc-ASAN compliant
