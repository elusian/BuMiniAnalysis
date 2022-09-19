# BuMiniAnalysis

## To produce the ntuple

Create a CMSSW area with the right version. Copy content of ntuwrite into src and compile.

After modifying the cfg with the right files and GT

```
cmsRun ntuwrite/cfg_MC_mini.py
```

Output (ntu.root by default) is produced in the current folder.

## To analyze the ntuple

Create a CMSSW area with the right version. Copy content of nturead into src and compile.
Create a file with all the paths of all the ntu files (from previous step) with an y in front
E.g.
```
$ cat ntu.list
y path/to/ntu1.root
y path/to/ntu2.root
y path/to/ntu3.root
```
Then run

```
pdTreeAnalyze ntu.list his.root -v outputFile second.root -v histoMode RECREATE
```
You can also add `-n nEvents` and `-s nSkip` to the command to limit the events you are running on

## To produce the plots

Create two folders called `plotSvt` and `plotPsi2S` (no, the macros won't do it).
You should have two different files from the ntuple analysis (second.root), place them in two folders called <basename>_base e <basename>_mkFit (basename was meant to be used when you have multiple comparison, e.g. noPU and PU)
Call
```
root -l -q -b drawSvtVars.C -- '"<basename>"'
root -l -q -b drawPsi2SVars.C -- '"<basename>"'
```
If your comparison is not mkFit vs no mkFit you probably need to change the labels inside the macro. Sorry for the mess
