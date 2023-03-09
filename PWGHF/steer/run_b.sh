o2-analysis-hf-task-b0 -b --configuration json://configuration.json |\
o2-analysis-hf-candidate-selector-b0-to-d-pi -b --configuration json://configuration.json |\
o2-analysis-hf-candidate-creator-b0 -b --configuration json://configuration.json |\
o2-analysis-hf-candidate-selector-dplus-to-pi-k-pi -b --configuration json://configuration.json |\
o2-analysis-hf-candidate-creator-3prong -b --configuration json://configuration.json |\
o2-analysis-hf-track-index-skim-creator -b --configuration json://configuration.json |\
o2-analysis-hf-track-to-collision-associator -b --configuration json://configuration.json |\
  o2-analysis-pid-tpc-full -b --configuration json://configuration.json |\
  o2-analysis-pid-tof-full -b --configuration json://configuration.json |\
  o2-analysis-pid-tpc-base -b --configuration json://configuration.json |\
  o2-analysis-pid-tof-base -b --configuration json://configuration.json |\
o2-analysis-trackselection -b --configuration json://configuration.json |\
o2-analysis-timestamp -b --configuration json://configuration.json |\
o2-analysis-track-propagation -b --configuration json://configuration.json |\
o2-analysis-collision-converter -b --configuration json://configuration.json --aod-file AO2D.root


