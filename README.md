# BCB420_GEO_GENE

## Usage
The usage is explained in `./R/main.R`. The main step is to create a JSON file as described below, and edit the parameters at the top of the R script.

Parameters to edit:
```
# json filename
inputFile <- "sample.json"

# output filename for RData
outputFile <- "sampletable.RData"

# log file
logFile <- "sample.log"
```

## JSON Format

Each experiment (comparison) is considered to be a JSON object, with the parameters "geoSeries", "description", "controlSamples", and "experimentalSamples". The values for each parameter should be a string. "geoSeries" describes the GEO Series identifier. "description" should be a unique description or descriptor for the experiment. "controlSamples" defines the normal sample(s), and "experimentalSamples" defines the experimental sample(s).

The JSON file is an array of experiments.

Sample JSON input:
```
[{
		"geoSeries": "GSE35330",
		"description": "Experiment 1",
		"controlSamples": ["GSM866121", "GSM866122", "GSM866124", "GSM866126"],
		"experimentalSamples": ["GSM866111", "GSM866116", "GSM866118", "GSM866119"]
	}
]
```
