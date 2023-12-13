o2-analysis-dq-table-maker -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-pid-tpc-full -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-pid-tof -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-pid-tof-full -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-pid-tof-beta -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-ft0-corrected-table -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-event-selection -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-trackselection -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-multiplicity-table -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-centrality-table -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-tracks-extra-converter -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-pid-tpc-base -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-pid-tof-base -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-timestamp -b --configuration json://configs/configTableMakerDataRun3.json |\
o2-analysis-track-propagation -b --configuration json://configs/configTableMakerDataRun3.json --aod-file @input_DataAOD_PbPb.txt --aod-writer-json configs/writerConfiguration_reducedEvent.json
