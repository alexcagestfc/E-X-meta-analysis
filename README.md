# E-X-meta-analysis
scripts developed for studying Engin-X instrument's beam line performance over time

EnginXMetaAnalysisBlocks  contains user defined functions that other programs draw from
### EXMAB will need to be edited to have correct file paths
file paths are expected at following points:

Section 1 script 1- identifyexperiments()
  script reads though a file called "SUMMARY.txt". You'll need to download file and copy file path to the script's 3rd line (line: 90)

Section 3 script 1- readexperiment(wkspce,exprimentaldict)
  script requres a file to save information to. This can be made by saving an empty notepad doccument and copying its file path with name to 2nd and 5th line in script (lines 292 and 295 in current version)
  note: change the file from ".txt" to ".csv" when entering file path

analysiscript when ran will call upon EXMAB to identify appropriate experiments, load experimental data, access loaded data, and analyse trends 

ranged analysis will assess performance for a predetermined range of energy levels, also using EXMAB
