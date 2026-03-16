# --- Permanent Project-Specific R Configuration ---

# 1. Force R to use wget (solves network download errors)
options(download.file.method = "wget")

# 2. Activate the renv library (solves R path conflicts)
# This will now succeed because .Renviron has fixed the rpm error.
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}
