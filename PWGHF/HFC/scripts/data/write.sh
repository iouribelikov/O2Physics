OPT='-b --configuration json://write.json --readers 1'
#OPT='-b --configuration json://write.json --readers 1 --aod-memory-rate-limit 471859200 --shm-segment-size 3758096384 --time-limit 12000'

o2-analysis-hf-pid-creator $OPT |\
o2-analysis-pid-tof-merge $OPT |\
o2-analysis-ft0-corrected-table $OPT |\
o2-analysis-tracks-extra-v002-converter $OPT |\
o2-analysis-multcenttable $OPT |\
o2-analysis-event-selection-service $OPT |\
o2-analysis-propagationservice $OPT |\
o2-analysis-pid-tpc-service $OPT |\
o2-analysis-hf-correlator-dplus-dplus-reduced $OPT |\
o2-analysis-hf-candidate-selector-dplus-to-pi-k-pi $OPT |\
o2-analysis-hf-candidate-creator-3prong $OPT --aod-file @input_data.txt --aod-writer-json OutputDirector.json --aod-parent-access-level 1 --aod-parent-base-path-replacement "alien:///alice/cern.ch/user/a/alihyperloop/jobs/0080/hy_809900/;"

