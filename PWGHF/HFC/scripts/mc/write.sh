OPT='-b --configuration json://write.json --readers 1'
#OPT='-b --configuration json://write.json --readers 1 --time-limit 600 --aod-memory-rate-limit 471859200 --shm-segment-size 3758096384'

o2-analysis-event-selection-service $OPT |\
o2-analysis-propagationservice $OPT |\
o2-analysis-tracks-extra-v002-converter $OPT |\
o2-analysis-multcenttable $OPT |\
\
o2-analysis-mccollision-converter $OPT |\
o2-analysis-ft0-corrected-table $OPT |\
\
o2-analysis-hf-pid-creator $OPT |\
o2-analysis-hf-mc-pid-tof $OPT |\
o2-analysis-pid-tpc-service $OPT |\
\
o2-analysis-hf-candidate-creator-3prong $OPT |\
o2-analysis-hf-candidate-selector-dplus-to-pi-k-pi $OPT |\
\
o2-analysis-hf-correlator-dplus-dplus-reduced $OPT \
 --aod-file @input_data.txt --aod-writer-json OutputDirector.json --aod-parent-access-level 1 --run
