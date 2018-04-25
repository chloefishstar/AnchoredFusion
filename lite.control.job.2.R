
## For job-cpu control
## To be used after for loop

        nruns = length(list.files('./', '_job_', full=F))
        while (nruns >0){
                Sys.sleep(30)
                nruns = length(list.files('./', '_job_', full=F))
        }

