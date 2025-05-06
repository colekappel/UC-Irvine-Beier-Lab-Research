# UC-Irvine-Beier-Lab-Research-LinearTrack

- **AnalysisFunctionPlaceCells.m** is the program that calculates stability, info score, sparsity, amplitude and event rate for just place cells - this program also excludes cells with low event rates/ spike data, **AnalysisFunction.m** is the same as **AnalysisFunctionPlaceCells.m** but it runs the function on all cells and does not exclude cells with low event rates/ spike data
- run **AnalysisProgramDay1_PlaceCells.m** to run **AnalysisFunctionPlaceCells.m** on all the mice - this will be the results for only place cells with exclusions of low event rate cells and run **AnalysisProgramDay1.m** to run **AnalysisFunction.m** on all mice - which will generate results for all cells in all mice without exclusions of cells with low event rates/ spike data
- after running either **AnalysisProgramDay1_PlaceCells.m** or **AnalysisProgramDay1.m**, run **EndBinsRemoved_AllPlotsAndStats.m** to get the statistics (lme models, t test, ks test, box plots, cdf plots, and heat maps of place cells) - the statistics compare the healthy mice to the mice in the 5xFAD+ group

- **infoScore.m, getEventRatePerBinMatx.m, randoSpikeMatrix.m and sparsityFunc.m** are just functions used by **AnalysisFunctionPlaceCells.m and AnalysisFunction.m**
