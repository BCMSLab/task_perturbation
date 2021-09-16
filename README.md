# task_perturbation
Code for the article "Disruption of cancer cell functions by task-specific drug perturbations"

To run this code,

- Build a docker image using the `Dockerfile` or pull the latest image [bcmslab/task_perturbation](https://hub.docker.com/r/bcmslab/task_perturbation)
- Run the three R scripts in `code/` in this order
    1. `code/prepart_data.R`: downloads the required data from source
    2. `code/analysis.R`: excutes the code to perform the analysis and save the output
    3. `code/output.R`: genrates the figures and graphs
