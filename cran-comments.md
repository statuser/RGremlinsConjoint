## Test environments
* local OS X install, R 4.1.2
* MacOS 10.15.7 (on GitHub Actions), R 4.1.0
* ubuntu 20.04 (on GitHub Actions), R 4.1.0
* ubuntu 20.04 (on GitHub Actions), R 2021-06-03 r80451
* windows-latest (on GitHub Actions), R 4.1.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

NOTES

Re-submission

* Fixed - Please rather use the Authors@R field and declare Maintainer, Authors and Contributors with their appropriate roles with person() calls.

* Fixed - Please do not start the description with "This package", package name, title or similar.

* Fixed - Please always write package names, software names and API (application programming interface) names in single quotes in title and description. e.g: --> 'gremlins'

* Fixed - If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file ...

* Fixed - You write information messages to the console that cannot be easily suppressed. (Added verbose option)



## Downstream dependencies
There are no downstream dependencies
