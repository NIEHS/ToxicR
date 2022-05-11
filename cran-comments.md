## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

checking installed package size ... NOTE
  installed size is  6.5Mb
  sub-directories of 1Mb or more:
    doc    1.4Mb
    libs   4.6Mb

The size of the the libs 4.6MB is based upon the compiler.  For macOS and Windows, it is usually around 5Mb. On Linux, it can be much larger ~ 100mb. 
This size can be reduced by non-standard flags (e.g., -Os), which are compiler specific. 
