# 

system("sed 's/: /=/g' config.txt | sed 's:{::' | sed 's:}::' | sed 's/, /; /g' > config.R")
source("config.R")

# Step 1...
# Step 2...

cat("==========\n", date(), ": Anchored-Fusion completed successfully.\n==========\n")
## END

