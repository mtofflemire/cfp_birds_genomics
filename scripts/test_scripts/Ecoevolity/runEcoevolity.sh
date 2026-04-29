#running the analysis
ecoevolity ecoevolity-config.yml #may get an error
ecoevolity --relax-missing-sites ecoevolity-config.yml




#summarize and plotting
ecoevolity --relax-triallelic-sites --relax-missing-sites ecoevolity-config.yml
pyco-sumchains -s 100 ecoevolity-config-state-run-1.log ecoevolity-config-state-run-2.log
pyco-sumchains -s 100 ecoevolity-config-state-run-1.log ecoevolity-config-state-run-2.log > pyco-sumchains-table.txt
sumcoevolity -b 100 -c ecoevolity-config.yml -n 1000000 ecoevolity-config-state-run-1.log ecoevolity-config-state-run-2.log





#plotting posterior probabilities
pyco-sumevents sumcoevolity-results-nevents.txt
pyco-sumtimes -b 100 -z ecoevolity-config-state-run-1.log ecoevolity-config-state-run-2.log



