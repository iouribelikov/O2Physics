#
#OPT='-b --configuration json://myconfig.json --readers 3 --shm-segment-size 3000000000 --aod-memory-rate-limit 471859200'
OPT='-b --configuration json://myconfig.json --readers 2'

o2-analysis-tracks-extra-v002-converter $OPT |\
o2-analysistutorial-task-yann $OPT --aod-file @input_data.txt

