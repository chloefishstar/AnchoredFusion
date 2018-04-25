
## for job-cpu control
## to be used in for loop

       system(paste("bash _job_.", i, " &", sep=''))

                Sys.sleep(0.1)
		if (!exists('nruns')) {nruns=1}
                while (nruns >= as.numeric(cpuMax) ){
                        Sys.sleep(30)
                        (nruns = length(list.files('./', '_job_', full=F)))
                }
